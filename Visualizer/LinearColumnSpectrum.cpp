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
	fExp = vz->genExpFreqLimits(numBands, vz->freqMin, vz->freqMax, vz->freq / 2);
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

void LinearColumnSpectrum::nextFrame()
{
	vz->source->getSample(intermediate);
	Utils::shortToFloat(readyOutBuffer, intermediate, bufferSize);

	vz->fft->winBlackman(readyOutBuffer, fftSize);

	dbLvl = vz->getDbLevel(readyOutBuffer, fftSize); //, dbMin, dbMax);
	hisDbLevel->push(dbLvl);
	silent = Utils::isSilence(hisDbLevel, vz->dbMin);

	//printf("silent: %d\t dbLvl: %f\n", (int)silent, dbLvl);

	vz->fft->process(fftDataBuffer, readyOutBuffer, fftSize);
	vz->scaleFft(fftDataBuffer, maxFftIndex);

	float avgBandEnergy = 0.0f;
	float avgHistoryEnergy = 0.0f;
	float old;
	
	for (unsigned int i = 0; i < numBands; i++) 
	{
		// calculate the average energy for the current band
		Freq f = fExp[i];

		int minFreqBandIndex = Utils::freqToIndex(vz->freq, f.min, maxFftIndex - 1);
		int maxFreqBandIndex = Utils::freqToIndex(vz->freq, f.max, maxFftIndex - 1);
		int delta = maxFreqBandIndex - minFreqBandIndex;

		avgBandEnergy = 0.0f;
		for (int e = minFreqBandIndex; e <= maxFreqBandIndex; e++)
			avgBandEnergy += fftDataBuffer[e];
		avgBandEnergy /= (delta + 1);

		if (silent)
			avgBandEnergy = 0;

		avgHistoryEnergy = hisAverage->avg();
		hisAverage->push(avgBandEnergy);

		// ---------------------------------

		int height = (int)round(winHi * Utils::norm(avgBandEnergy, 0, 50));
		int step = (int)Utils::lerp(prevBandEnergies[i], height, 0.3f);
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