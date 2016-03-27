#include "LinearColumnSpectrum.h"


LinearColumnSpectrum::LinearColumnSpectrum()
{ }

LinearColumnSpectrum::~LinearColumnSpectrum()
{ }

void LinearColumnSpectrum::init(const Visualizer* vzInstance)
{
	vz = (Visualizer*)vzInstance;

	window = new sf::RenderWindow(sf::VideoMode(850, 450), "wub wub wuuub-b-b");
	backgroundColor = sf::Color::Black;
	
	bufferSize = vz->fft->getBufferSize();
	fftSize = vz->fft->getFFTSize();
	maxFftIndex = vz->fft->getMaxFftIndex();

	intermediate = new unsigned char[88200];
	readyOutBuffer = new float[fftSize];
	fftDataBuffer = new float[fftSize];

	for (unsigned int i = 0; i < fftSize; i++)
	{
		readyOutBuffer[i] = 0.0f;
		fftDataBuffer[i] = 0.0f;
	}

	maxHistoryEntries = 30;
	numBands = 32;
	silent = false;

	prevBandEnergies = new float[numBands];
	for (unsigned int i = 0; i < numBands; i++)
		prevBandEnergies[i] = 0.0f;

	hisDbLevel = new History(maxHistoryEntries);
	hisAverage = new History(20);

	winWi = window->getSize().x;
	winHi = window->getSize().y;
	rectWi = (float)(winWi / numBands);
	sf::Color barColor = sf::Color(202, 44, 63);

	for (unsigned int i = 0; i < numBands; i++)
	{
		sf::RectangleShape *rect = new sf::RectangleShape();
		rect->setSize(sf::Vector2f(rectWi, 20));
		rect->setPosition(sf::Vector2f(i * rectWi + (2 * i), winHi - rect->getSize().y));
		rect->setFillColor(barColor);
		rects.push_back(rect);
	}

	fLin = vz->genLinFreqLimits(numBands, vz->freqMin, vz->freqMax);
	fLog = vz->genLogFreqLimits(numBands, vz->freq);
	fExp = vz->genExpFreqLimits(numBands, vz->freqMin, vz->freqMax, vz->freq / 2);

	lp = new Biquad();
	bp = new Biquad();
	pf = new Biquad();

	float ff = 400.0f;
	lp->setBiquad(bq_type_lowpass, (float)(ff / (float)(vz->freq)), 0.707, 0);
	bp->setBiquad(bq_type_bandpass, (float)(ff / (float)(vz->freq)), 0.9, 0);
	pf->setBiquad(bq_type_peak, (float)(ff / (float)(vz->freq)), 0.7071, 3);

	// beat detection initialization

	maxBandEntries = 43;
	beatCount = 0;
	beatClock.restart();

	hisEnergy = new History(maxBandEntries);
	
	// individual band histories
	for (unsigned int b = 0; b < numBands; b++)
	{
		bandHistories.push_back(new History(maxBandEntries));
	}

	// random thing
}

void LinearColumnSpectrum::applyFilter(float *inData, size_t inLen)
{ 
	for (unsigned int i = 0; i < inLen; i++)
	{
		lp->process(inData[i]);
		//bp->process(inData[i]);
		//pf->process(inData[i]);
	}
}

void LinearColumnSpectrum::start()
{
	vz->source->start();

	while (window->isOpen())
	{
		sf::Event event;
		while (window->pollEvent(event))
		{
			if (event.type == sf::Event::Closed)
				window->close();

			//if (event.type == sf::Event::MouseButtonReleased) { }
		}

		frameDelta = deltaClock.restart();
		deltaTime = frameDelta.asSeconds();

		nextFrame();

		// do the drawing
		window->clear(backgroundColor);
		
		for (auto* r : rects)
			window->draw(*r);
		window->draw(dbMeter);

		window->display();

		Sleep(vz->sleepTime);
	}

	vz->source->stop();
	stop();
}

float calc_c(float variance) //, float f)
{
	// TODO: add comments so this is not just a bunch of magic numbers

	return ((-1.0f) * 0.0025714f * variance) + 1.3142857f;
	//return ((-1.0f) * 0.0000002f * variance) + 1.3142857f;
	//return ((-1.0f) * 0.0025714f * variance) + 1.5142857f;
	//return ((-1.0f) * 0.0025714f * variance) + 5.559142857f;
}

void LinearColumnSpectrum::nextFrame()
{
	vz->source->getSample(intermediate);
	Utils::shortToFloat(readyOutBuffer, intermediate, bufferSize);

	vz->fft->winBlackman(readyOutBuffer, fftSize);

	dbLvl = vz->getDbLevel(readyOutBuffer, fftSize); //, dbMin, dbMax);
	hisDbLevel->push(dbLvl);
	silent = Utils::isSilence(hisDbLevel, vz->dbMin);
	//printf("silent: %d\t dbLvl: %f\n", (int)silent, dbLvl);

	//applyFilter(readyOutBuffer, fftSize);

	vz->fft->process(fftDataBuffer, readyOutBuffer, fftSize);

	applyFilter(fftDataBuffer, maxFftIndex);

	// do some beat detection woopwoop

	// maybe change out maxFftIndex for fftSize?
	float energy = 0.0f;
	for (unsigned int k = 0; k < maxFftIndex; k++)
		energy += fftDataBuffer[k];
	//energy /= maxFftIndex;

	float avgEnergy = hisEnergy->avg();
	float varEnergy = hisEnergy->var(avgEnergy);
	float C = calc_c(varEnergy);
	hisEnergy->push(energy);

	float shit = C * avgEnergy;
	if (energy > shit) 
	{
		if (beatClock.getElapsedTime().asMilliseconds() > 150)
		{
			beatClock.restart();
			beatCount++;
			printf("BEAT %d\n", beatCount);
			// random things
			// usage: dist(md)
			std::random_device rd;
			std::mt19937 mt(rd());
			std::uniform_int_distribution<int> dist(0, 128);
			backgroundColor = sf::Color(dist(mt), dist(mt), dist(mt)); //getRandomColor(dist, mt);
		}
	}

	float avgBandEnergy = 0.0f;
	float avgHistoryEnergy = 0.0f;
	float old;

	// do some drawing with the individual bands
	vz->scaleFft(fftDataBuffer, maxFftIndex);

	for (unsigned int i = 0; i < numBands; i++) 
	{
		// calculate the average energy for the current band
		Freq f = fLin[i];

		int minFreqBandIndex = Utils::freqToIndex(vz->freq, f.min, maxFftIndex - 1);
		int maxFreqBandIndex = Utils::freqToIndex(vz->freq, f.max, maxFftIndex - 1);
		int delta = maxFreqBandIndex - minFreqBandIndex;

		avgBandEnergy = 0.0f;
		for (int e = minFreqBandIndex; e <= maxFreqBandIndex; e++)
			avgBandEnergy += fftDataBuffer[e];
		avgBandEnergy /= (delta + 1);

		if (silent)
			avgBandEnergy = 0;

		bandHistories[i]->push(avgBandEnergy);
		float avgBandHistoryEnergy = bandHistories[i]->avg();

		// ---------------------------------

		int height = (int)round(winHi * Utils::norm(avgBandEnergy, 0, 100));
		int step = 0;
		float pe = prevBandEnergies[i];
		if (pe < height)
			step = (int)Utils::lerp(prevBandEnergies[i], height, 0.85f);
		else
			step = (int)Utils::lerp(prevBandEnergies[i], height, 0.5f);
		prevBandEnergies[i] = step;

		height = max(min(step, winHi), 1);
		//height = max(min(height, winHi), 1);
				
		old = height;
		height = (height + prevRectHeight) / (i == 0 ? 1 : 2);
		prevRectHeight = old;

		currPos.x = i * rectWi + (2 * i);
		currPos.y = winHi - height;
		currSize.y = height;
		currSize.x = rectWi;

		rects[i]->setPosition(currPos);
		rects[i]->setSize(currSize);
	}
}

void LinearColumnSpectrum::stop()
{
	if (window->isOpen())
	{
		window->close();
	}
	// also delete other resources here
}