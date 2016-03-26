#include "main.h"

#include "fftw3.h"
#include "Biquad.h"
#include "al.h"
#include "alc.h" 
#include "SerialClass.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <list>
#include <vector>
#include <random>

#include <Windows.h>

#include <SFML/Graphics.hpp>
#include <SFML/Audio.hpp>

//struct Freq
//{
//	float min;
//	float max;
//};
// typedef std::vector<freq> FreqBand;

void check_al_error()
{
	ALCenum err = alGetError();
	if (err != AL_NO_ERROR)
	{
		printf("Error code: %d", err);
		system("PAUSE");
	}
}

float norm(float x, float min, float max)
{
	return (x - min) / (max - min);
}

float lerp(float start, float end, float t) 
{
	return (1.0f - t) * start + t * end;
}

int freq_to_index(size_t sample_rate, float freq, size_t max_fft_idx)
{
	float ny = sample_rate / 2.0f;
	return (int)roundf((freq / ny) * max_fft_idx);
}

std::vector<Freq> gen_freq_limits(int n, int min, int max, size_t sample_rate)
{
	const int delta = max - min;
	const int freqDelta = (int)(delta / n);
	std::vector<Freq> freqs;

	for (unsigned int i = 0; i < n; i++)
	{
		Freq f;

		f.min = floorf(min + (i * freqDelta)); // freq_to_index(sample_rate, 
		f.max = floorf(min + ((i + 1) * freqDelta)); // freq_to_index(sample_rate, 

		freqs.push_back(f);
	}

	return freqs;
}

std::vector<Freq> gen_log_freq_limits(int n, size_t sample_rate)
{
	std::vector<Freq> freqs;

	for (unsigned int i = 0; i < n; i++)
	{
		float avg = 0.0f;
		int lowFreq = 0;

		if (i == 0)
			lowFreq = 0;
		else
			lowFreq = (int)((sample_rate / 2) / (float)pow(2, n - i));

		int hiFreq = (int)((sample_rate / 2) / (float)pow(2, (n - 1) - i));

		int lowBound = lowFreq; //freq_to_index(sample_rate, lowFreq, n);
		int highBound = hiFreq; //freq_to_index(sample_rate, hiFreq, n);

		Freq f;

		f.min = lowBound;
		f.max = highBound;

		freqs.push_back(f);
	}

	return freqs;
}

std::vector<Freq> gen_exp_freq_limits(int n, int min, int max, size_t sample_rate)
{
	std::vector<Freq> freqs;

	auto fx = [](int x, int n, int min, int max)
	{
		//double ex = 1.05;
		//return ((max - min) / pow(ex, n)) * pow(ex, x) + min;
		double ex = 1.2f; //2.718281828459; //1.5; //sqrt(2);
		return (float)((max - min) / pow(n, ex)) * pow(x, ex) + min;
	};

	for (unsigned int i = 0; i < n; i++)
	{
		Freq f;

		f.min = ceilf(fx(i, n, min, max)); //freq_to_index(sample_rate, floorf(fx(i, n, min, max)), n);
		f.max = ceilf(fx(i + 1, n, min, max)); //freq_to_index(sample_rate, floorf(fx(i + 1, n, min, max)), n);

		freqs.push_back(f);
	}

	return freqs;
}

float calc_channel_energy(float *in, size_t len)
{
	float res = 0.0f;

	for (unsigned int i = 0; i < len; i++)
	{
		res += in[i] * in[i];
	}

	return res;
}

float calc_c(float variance, float f)
{
	return ((-1.0f) * 0.0025714f * variance) + f; //1.3142857f;
	//return ((-1.0f) * 0.0000002f * variance) + f; //1.3142857f;
	//return ((-1.0f) * 0.0025714f * variance) + 1.5142857f;
	//return ((-1.0f) * 0.0025714f * variance) + 5.559142857f;
}

void blackman_window(float *in, size_t sampleLen)
{
	float a0 = (1 - 0.16f) / 2;
	float a1 = 0.5f;
	float a2 = 0.16f / 2;

	for (unsigned int i = 0; i < sampleLen; i++)
	{
		float mult = a0 - a1*cosf((2.0f * M_PI * i) / (sampleLen - 1))
			+ a2*cosf((4.0f * M_PI * i) / (sampleLen - 1));
		in[i] *= mult;
	}
}

void hann_window(float *in, size_t sampleLen)
{
	for (unsigned int i = 0; i < sampleLen; i++)
	{
		float mult = 0.5f * (1.0f - cosf(2.0f * M_PI * i / (sampleLen - 1)));
		in[i] *= mult;
	}
}

void hamming_window(float *in, size_t sampleLen)
{
	float a = 0.54f;
	float b = 0.46f;

	for (unsigned int i = 0; i < sampleLen; i++)
	{
		float mult = a - (b * cosf(2.0f * M_PI * i / (sampleLen - 1)));
		in[i] *= mult;
	}
}

void harris_window(float *in, size_t sampleLen)
{
	float a0 = 0.35875f;
	float a1 = 0.48829f;
	float a2 = 0.14128f;
	float a3 = 0.01168f;

	for (unsigned int i = 0; i < sampleLen; i++)
	{
		float mult = a0 
			- (a1 * cosf(2.0f * M_PI * i / (sampleLen - 1))) 
			+ (a2 * cosf(4.0f * M_PI * i / (sampleLen - 1)))
			+ (a3 * cosf(6.0f * M_PI * i / (sampleLen - 1)));
		in[i] *= mult;
	}
}

void short_to_float(float *out, unsigned char *in, size_t len)
{
	for (unsigned int i = 0; i < len; i++)
	{
		int n = (in[i] << 8 | in[i + 1]);
		out[i / 2] = (float)n / 32768.0f; // signed short max value
	}
}

bool is_silence(std::list<float> *his, int min_db)
{
	bool ret = true;

	for (std::list<float>::iterator it = his->begin(); it != his->end(); it++)
	{
		if ((int)round(*it) > min_db)
		{
			ret = false;
		}
	}

	return ret;
}

void push_history(float val, std::list<float> *his, size_t maxlen)
{
	if (his->size() == maxlen)
		his->pop_back();

	his->push_front(val);
}

float sum_history(std::list<float> *his)
{
	float res = 0.0f;

	for (std::list<float>::iterator it = his->begin(); it != his->end(); it++)
	{
		res += *it;
	}

	return res;
}

float avg_history(std::list<float> *his)
{
	float res = 0.0f;
	
	for (std::list<float>::iterator it = his->begin(); it != his->end(); it++)
	{
		res += *it;
	}

	res /= his->size();
	return res;
}



float average(float in[], size_t len)
{
	float ret = .0f;
	for (int i = 0; i < len; i++) {
		ret += in[i];
	}
	return ret /= len;
}

float average2(float in[], size_t len)
{
	float ret = .0f;
	for (int i = 0; i < len; i++) {
		ret += in[i];
	}
	if (len > 0) {
		ret /= len;
	}
	return ret;
}

float average2(std::list<float> *in)
{
	float ret = .0f;
	if (in->size() == 0) {
		return ret;
	}
	for (std::list<float>::iterator it = in->begin(); it != in->end(); it++) {
		ret += *it;
	}
	return ret /= in->size();
}

float variance(float in[], size_t len, float energy)
{
	float ret = .0f;
	for (int i = 0; i < len; i++) {
		ret += powf(in[i] - energy, 2);
	}
	return ret /= len;
}

float variance(std::list<float> *in, float energy)
{
	float ret = .0f;
	for (std::list<float>::iterator it = in->begin(); it != in->end(); it++) {
		ret += powf(*it - energy, 2);
	}
	return ret /= in->size();
}

// why do i even have this

/*
TODO:
	- Use the actual frequency of the recording device
	- Create a GetBufferData interface so I can use the winapi more easily
	- Split up the rendering into another class or interface for simpler code
	- Put all the math into one file
	- just clean up the code dammit
	- Add proper error handling for OpenAL

	gamedev.net/page/resources/_/technical/math-and-physics/beat-detection-algorithms-r1952
*/

int main()
{
	Visualizer vz(FFTSize::FFT1024);
	vz.init();

	return 0;
}

int main2()
{
	const unsigned int bufferSize = 2048; // 2048; //2048;
	// 2 8-bit samples -> 1 float
	// which means the length is half the size
	const unsigned int fftSize = bufferSize / 2; // 1024; //1024;
	const unsigned int numBands = 30; // 12;
	const unsigned int freq = 44100;
	const long actualSleepTime = 40; // ms

	int min_db = 0;
	int max_db = -90;
	int min_freq = 20;
	int max_freq = 4000;

	unsigned int maxFftIndex = (fftSize / 2) + 1;

	float fftBuffer[fftSize];
	float mags[fftSize];
	//bool recording = true;
	//ALCint prev = 0;
	ALCint dataLen = 0;
	ALCenum err;
	unsigned char buffer[88200];
	unsigned int buffersPerSec = 20; //round(freq / bufferSize);

	bool useSerial = true;
	char* port = "COM3"; /* \\\\.\\ */
	char serialBuffer[31 * 3];
	size_t serialBufferLen = 31 * 3;
	Serial* sc = new Serial(port);

	if (!sc->IsConnected()) {
		printf("Could not connect to serial port. LED disabled");
		useSerial = false;
	}
	else {
		printf("Serial connection opened successfully.\n");
	}

	fftwf_plan plan;
	fftwf_complex* in;
	fftwf_complex* out;
	
	std::vector<Freq> fLin = gen_freq_limits(numBands, min_freq, max_freq, freq / 2);
	std::vector<Freq> fLog = gen_log_freq_limits(numBands, freq);
	std::vector<Freq> fExp = gen_exp_freq_limits(numBands, min_freq, max_freq, freq / 2);

	// www.fftw.org/doc/One_002dDimensional-DFTs-of-Real-Data.html#One_002dDimensional-DFTs-of-Real-Data
	//in = (fftwf_complex*)fftw_alloc_complex(sizeof(fftwf_complex) * fftSize);
	out = (fftwf_complex*)fftw_alloc_complex(sizeof(fftwf_complex) * fftSize);
	plan = fftwf_plan_dft_r2c_1d(fftSize, fftBuffer, out, FFTW_ESTIMATE);

	const ALchar *pDeviceDefault = alcGetString(
		NULL, 
		ALC_CAPTURE_DEVICE_SPECIFIER); // get default capture device

	if ((err = alGetError()) != AL_NO_ERROR)
	{
		printf("Error code: %d", err);
	}

	// blitzbasic.com/Community/posts.php?topic=90830
	ALCdevice *pCaptureDevice = alcCaptureOpenDevice(
		pDeviceDefault,
		freq,               // samples per second
		AL_FORMAT_MONO16,  // one channel with two 8-bit samples
		bufferSize); // sample_rate * res * num_tracks * sec

	if ((err = alGetError()) != AL_NO_ERROR)
	{
		printf("Error code: %d", err);
	}

	if (!pCaptureDevice)
	{
		printf("Could not create capture device.\n");
		exit(1);
	}

	printf("Using default recording device:\n\t%s\n\n", pDeviceDefault);
	printf("Sleep time: %d ms\n", buffersPerSec);
	printf("Operating at frequency: %d\n", freq);
	printf("Format: AL_FORMAT_MONO16\n");
	printf("Capture buffer size: %d\n", bufferSize);
	printf("FFT buffer size: %d\n", fftSize);

	printf("\nNumber of bands: %d\n", numBands);
	printf("Min db: %d\n", min_db);
	printf("Max db: %d\n", max_db);
	printf("Min freq: %d\n", min_freq);
	printf("Max freq: %d\n", max_freq);

	printf("\nStarting capture device...\n");
	alcCaptureStart(pCaptureDevice);
	
	if ((err = alGetError()) != AL_NO_ERROR)
	{
		printf("Error code: %d", err);
	}

	// random things
	// usage: dist(md)
	std::random_device rd;
	std::mt19937 mt(rd());
	std::uniform_int_distribution<int> dist(0, 128);

	printf("Initalizing render window...\n");
	sf::RenderWindow window(sf::VideoMode(850, 450), "wub wub wuuub-b-b"); //550, 270
	std::vector<sf::RectangleShape> rects;
	sf::Vector2f vecPos;
	sf::Vector2f vecSize;
	vecSize.x = window.getSize().x / numBands;

	int win_hi = window.getSize().y;
	int win_wi = window.getSize().x;
	float rectWidth = 0.0f;
	
	for (unsigned int k = 0; k < numBands; k++)
	{
		rectWidth = window.getSize().x / numBands;
		float botPos = window.getSize().y;
		sf::RectangleShape rect;
		
		rect.setSize(sf::Vector2f(rectWidth, 20));
		rect.setPosition(sf::Vector2f(k * rectWidth + (2 * k), botPos - rect.getSize().y));
		rect.setFillColor(sf::Color(202, 44, 63));

		rects.push_back(rect);
	}

	sf::RectangleShape beatRect;
	beatRect.setSize(sf::Vector2f(10, 40));
	beatRect.setPosition(sf::Vector2f(0, 0));
	beatRect.setFillColor(sf::Color::Green);
	float beatHi = 5.0f;
	float prev_hi = 0.0f;

	std::list<float> eBandHistory[numBands];
	std::list<float> eDiffHistory[numBands];
	size_t maxHistorySamples = 1; //buffersPerSec; //42; //sleep; //43;

	//float lastAverage = 0.0f;
	//float lastEnergy = 0.0f;
	float C = 200.0f; //1.1f;//1.1f;
	bool isBeat = false;
	//float maxEnergy = 0;
	unsigned int beatCount = 0;
	
	float prevBandThing[numBands]; // used to band animation interpolation
	
	std::list<float> hisDbLevel;
	int hisDbLevelCount = buffersPerSec * 3; // should be aprox 1-3~ seconds
	float currDbLevel = 0; // the current aprox. db level of the current signal
	float avgDbLevel = 0;
	float prevAvgDbLevel = 0;
	bool isSilence = false;

	sf::Clock clock;
	sf::Time frameDelta;
	float dt = 0.0f;

	sf::Clock beatClock;
	sf::Time beatTime;
	float timeSinceLastBeat = 0.0f;

	// Beat detection
	std::list<float> beatHistory;
	int beatHistoryCount = buffersPerSec;
	float instantBeatEnergy = 0;
	float beatHistoryAverage = 0;
	int beatCount2 = 0;
	float beatWaitTime = 300; // ms

	float beatEnergies[numBands] = { 0.0f };
	bool beatMap[numBands] = { false };
	sf::Clock beatClocks[numBands] = { sf::Clock() }; // is this legal?

	beatClocks[4].restart();

	Biquad *lp = new Biquad();
	Biquad *bp = new Biquad();
	Biquad *pf = new Biquad();

	float ff = 250.0f;
	lp->setBiquad(bq_type_lowpass, (float)(ff / (float)(freq)), 0.707, 0);
	bp->setBiquad(bq_type_bandpass, (float)(ff / (float)(freq)), 0.9, 0);
	pf->setBiquad(bq_type_peak, (float)(ff / (float)(freq)), 0.7071, 3);

	beatClock.restart();

	float shitFuck = 1.3f; //142857f;
	float shitStep = 0.1f;
	sf::Color backgroundColor = sf::Color::Black;

	for (int i = 0; i < numBands; i++) {
		beatClocks[i].restart();
	}

	// real time streaming of the data or something
	// maybe just restart the whole thing because I
	// have no idea what I am doing at this point

	printf("Starting render loop...\n");
	while (window.isOpen())
	{
		//
		// SFML event polling
		//
		sf::Event event;
		while (window.pollEvent(event))
		{
			if (event.type == sf::Event::Closed)
				window.close();
			
			//if (event.type == sf::Event::MouseButtonReleased) {
			//}
		}
		//if (sf::Mouse::isButtonPressed(sf::Mouse::Button::Left)) shitFuck += shitStep;
		//if (sf::Mouse::isButtonPressed(sf::Mouse::Button::Right)) shitFuck -= shitStep;
		window.setTitle(std::to_string(shitFuck));

		//
		// get the delta time of the frame
		//
		frameDelta = clock.restart();
		dt = frameDelta.asSeconds();

		//
		// Get the samples
		// 
		alcGetIntegerv(
			pCaptureDevice, 
			ALC_CAPTURE_SAMPLES,
			(ALCsizei)sizeof(ALCint),
			&dataLen);

		if ((err = alGetError()) != AL_NO_ERROR)
		{
			printf("Error code: %d", err);
		}

		// if we don't have enough samples yet,
		// we just wait another frame and just 
		// proceed with the previous one
		if (dataLen >= bufferSize / 1)
		{
			alcCaptureSamples(
				pCaptureDevice,
				(ALCvoid*)buffer,
				bufferSize);

			if ((err = alGetError()) != AL_NO_ERROR)
			{
				printf("Error code: %d", err);
			}
		}
		else
		{
			// printf("Warning: dataLen < captureBufferLen\n");
		}

		//
		// convert the mono samples to float
		// 16-bit audio = two 8-bit samples
		//
		short_to_float(fftBuffer, (unsigned char*)buffer, bufferSize);
		
		//
		// Apply a window function
		//
		//hann_window(fftBuffer, fftSize);
		blackman_window(fftBuffer, fftSize);
		//hamming_window(fftBuffer, fftSize);
		//harris_window(fftBuffer, fftSize);

		//
		// Get the approx. dB level of the current signal
		//
		currDbLevel = 0;
		float db = 0.0f;
		for (unsigned int i = 0; i < fftSize; i++)
		{
			db = 20.0f * log10(fftBuffer[i]);
			currDbLevel += db;
		}
		currDbLevel /= fftSize;

		// clamp db level to range 0 to 90
		currDbLevel = max(min(currDbLevel, min_db), max_db) - max_db;

		push_history(currDbLevel, &hisDbLevel, hisDbLevelCount);
		isSilence = is_silence(&hisDbLevel, 0);
		if (hisDbLevel.size() < hisDbLevelCount)
		{
			isSilence = false;
		}
		
		for (unsigned int i = 0; i < fftSize; i++)
		{
			fftBuffer[i] = 
				//pf->process(
				//lp->process(
				//bp->process(
					fftBuffer[i]
					;
		}
		
		//
		// Execute the FFT plan
		//
		fftwf_execute(plan);

		//
		// Calculate the magnitude of the signal and store it in an array
		//
		float f = 0.0f;
		for (unsigned int i = 0; i < maxFftIndex; i++) 
		{
			f = sqrtf((out[i][0] * out[i][0]) + (out[i][1] * out[i][1]));
			//f = lp->process(f);
			mags[i] = f;
		}

		//
		// Perform scaling, either logarithmically or exponentially
		//
		for (unsigned int i = 0; i < maxFftIndex; i++) 
		{
			//mags[i] *= 20;
			//mags[i] = 20.0f * log10(mags[i]);
			mags[i] = sqrt(mags[i]) * 2;
			//mags[i] = sqrtf(mags[i]);
			//mags[i] = 10.0f * log10(mags[i]);
		}

		float prev = 0;
		float avgBandDb = 0;
		float avgBandEnergy = 0;
		float currBandHistoryAvg = 0;
		float currBandHistoryVar = 0;

		//
		//
		//
		int minFreqBandIndex = freq_to_index(freq, min_freq, maxFftIndex - 1);
		int maxFreqBandIndex = freq_to_index(freq, max_freq, maxFftIndex - 1);
		float sum = 0;
		for (int e = minFreqBandIndex; e <= maxFreqBandIndex; e++)
		{
			sum += mags[e] * mags[e];
		}
		sum /= maxFftIndex;

		beatHistoryAverage = avg_history(&beatHistory);
		beatHistoryAverage = max(beatHistoryAverage, 0);

		//sum /= maxFftIndex;
		push_history(sum, &beatHistory, beatHistoryCount);

		float beatVariance = 0;
		for (float &s : beatHistory) beatVariance += pow(s - beatHistoryAverage, 2);
		beatVariance /= beatHistory.size();
		float bC = 250.0f; //calc_c(beatEnergyAverage);
		

		sf::Int32 elap = beatClock.getElapsedTime().asMilliseconds();
		if (sum > (bC * beatVariance) && !isSilence && elap >= 300)
		{
			printf("BEAT! %d\n", beatCount2);
			beatCount2++;

			beatClock.restart();
		}

		for (unsigned int i = 0; i < numBands; i++) {
			// calculate the average energy for the current band
			Freq f = fLin[i];
			avgBandEnergy = 0.0f;

			int minFreqBandIndex = freq_to_index(freq, f.min, maxFftIndex - 1);
			int maxFreqBandIndex = freq_to_index(freq, f.max, maxFftIndex - 1);
			int delta = maxFreqBandIndex - minFreqBandIndex;

			for (int e = minFreqBandIndex; e <= maxFreqBandIndex; e++) {
				avgBandEnergy += mags[e];
			}
			avgBandEnergy /= (delta + 1);

			if (isSilence) {
				avgBandEnergy = 0;
			}

			float eHistoryAverage = avg_history(&eBandHistory[i]);

			sf::Int32 elap = beatClocks[i].getElapsedTime().asMilliseconds();
			//if (elap <= beatWaitTime) {
			//	beatMap[i] = false;
				//rects[i].setFillColor(sf::Color(202, 44, 63));
			//}
			//else if (eDiff2 > 0) {
			//	beatMap[i] = true;
				//rects[i].setFillColor(sf::Color::Magenta);
			//	beatClocks[i].restart();
			//}
			//else {
			//	beatMap[i] = false;
				//rects[i].setFillColor(sf::Color(202, 44, 63));
			//}

			push_history(avgBandEnergy, &eBandHistory[i], beatHistoryCount);
			//push_history(eConstantDiff, &eDiffHistory[i], beatHistoryCount);

			beatEnergies[i] = avgBandEnergy;
			
			// --------------------------------------------------------------

			int height = (int)round((win_hi / 1) * norm(avgBandEnergy, 0, 50)); //norm(avgBandDb, min_db, max_db));
			int step = (int)lerp(prevBandThing[i], height, 0.7f);
			prevBandThing[i] = step;

			// clamp to window heigth
			height = max(min(height, win_hi), 1);

			// simple smoothing using last rect as reference
			//float old = height;
			//height = (height + prev) / (i == 0 ? 1 : 2);
			//prev = old;

			// set the rectangle position and size
			vecPos.x = i * rectWidth + (2 * i);
			vecPos.y = win_hi - height;
			vecSize.y = height;

			rects[i].setPosition(vecPos);
			rects[i].setSize(vecSize);

		}
		
		// --------------------------------------------------------------

		int count = 0;
		int rangeStart = 0;
		int rangeEnd = 10;
		int thresh = 3;

		for (int j = rangeStart; j < rangeEnd; j++) {
			if (beatMap[j]) {
				count++;
			}
		}

		for (int j = rangeStart; j < rangeEnd; j++) {
			if (count >= thresh) {
				//rects[j].setFillColor(sf::Color::Magenta);
				// backgroundColor = sf::Color(dist(mt), dist(mt), dist(mt)); //getRandomColor(dist, mt);
				// printf("awjkdajwkld %d\n", j);
			}
			else {
				//rects[j].setFillColor(sf::Color(202, 44, 63));
			}
		}

		// AKKKKKK
		if (count >= thresh) {
			beatCount++;
			//backgroundColor = sf::Color(dist(mt), dist(mt), dist(mt)); //getRandomColor(dist, mt);
			printf("awjkdajwkld %d\n", beatCount);

		}

		if (useSerial) {

			//for (int k = 0; k < serialBufferLen; k += 3) {
			serialBuffer[0] = backgroundColor.r;
			serialBuffer[0 + 1] = backgroundColor.g;
			serialBuffer[0 + 2] = backgroundColor.b;
			//}

			sc->WriteData(serialBuffer, 3); //serialBufferLen);
		}

		// --------------------------------------------------------------

		for (unsigned int i = 0; i < numBands; i++) {

			float e = beatEnergies[i];
		}

		//
		// Do beat detection and prepare graphics for drawing
		//
		for (unsigned int i = 0; i < numBands; i++) 
		{
			/*
			// calculate the average energy for the current band
			Freq f = fExp[i]; // fLog, fLin
			avgBandEnergy = 0.0f;

			int minFreqBandIndex = freq_to_index(freq, f.min, maxFftIndex - 1);
			int maxFreqBandIndex = freq_to_index(freq, f.max, maxFftIndex - 1);
			int delta = maxFreqBandIndex - minFreqBandIndex;

			for (int e = minFreqBandIndex; e <= maxFreqBandIndex; e++) 
			{
				avgBandEnergy += mags[e];
			}
			avgBandEnergy /= (delta + 1);

			// get the current history average for this band
			currBandHistoryAvg = avg_history(&eBandHistory[i]);

			// push the new avg onto into history
			push_history(avgBandEnergy, &eBandHistory[i], maxHistorySamples);
			
			// calculate the avg variance of the current subband's history
			currBandHistoryVar = 0;
			for (auto &s : eBandHistory[i]) currBandHistoryVar += pow(s - currBandHistoryAvg, 2);
			currBandHistoryVar /= eBandHistory[i].size();
			
			C = calc_c(currBandHistoryVar);
			if (i == 3)
			{
				window.setTitle("wubuwbubw - Variance: " 
					+ std::to_string(currBandHistoryVar) 
					+ " - C: " + std::to_string(C)
					+ " - Avg: " + std::to_string(currBandHistoryAvg));
			}

			if (avgBandEnergy > C * currBandHistoryAvg && !isSilence || false) 
			{
				//printf("Beat - %ld - i == %d\n", beats, i);
				//rects[i].setFillColor(sf::Color::Magenta);
				//beats++;
				//isBeat = true;
				
				int minKickIndex = freq_to_index(freq, 20, maxFftIndex - 1);
				int maxKickIndex = freq_to_index(freq, 100, maxFftIndex - 1);
				if (minFreqBandIndex >= minKickIndex 
					&& maxFreqBandIndex <= maxKickIndex 
					&& !isBeat) 
				{

					sf::Int32 elap = beatClock.getElapsedTime().asMilliseconds();
					if (elap >= 300) {
						//printf("Sub bass - - i == %d\n", i);
						rects[i].setFillColor(sf::Color::Magenta);
						beatCount++;
						isBeat = true;

						beatTime = beatClock.restart();
					}

				}

				minKickIndex = freq_to_index(freq, ff - 10, maxFftIndex - 1);
				maxKickIndex = freq_to_index(freq, ff + 10, maxFftIndex - 1);
				if (i >= minKickIndex && i <= maxKickIndex && !isBeat) 
				{
					//printf("Beat - %ld - i == %d\n", beatCount, i);
					//rects[i].setFillColor(sf::Color::Magenta);
					beatCount++;
					//isBeat = true;
				}
				minKickIndex = freq_to_index(freq, 220, maxFftIndex - 1);
				maxKickIndex = freq_to_index(freq, 240, maxFftIndex - 1);
				if (i >= minKickIndex && i <= maxKickIndex && !isBeat) 
				{
					// printf("Snare - %ld - i == %d\n", beatCount, i);
					//rects[i].setFillColor(sf::Color::Magenta);
					beatCount++;
					//isBeat = true;
				}
			}
			else 
			{
				rects[i].setFillColor(sf::Color(202, 44, 63));
				isBeat = false;
			}

			if (isSilence) avgBandEnergy = 0;
			

			// normalize the height of the rectangle to fit in the window
			int height = (int)round((win_hi / 1) * norm(avgBandEnergy, 0, 90)); //norm(avgBandDb, min_db, max_db));
			int step = (int)lerp(prevBandThing[i], height, 0.5f);
			prevBandThing[i] = step;

			// clamp to window heigth
			height = max(min(step, win_hi), 1);

			// simple smoothing using last rect as reference
			float old = height;
			height = (height + prev) / (i == 0 ? 1 : 2);
			prev = old;

			// set the rectangle position and size
			vecPos.x = i * rectWidth + (2 * i);
			vecPos.y = win_hi - height;
			vecSize.y = height;

			rects[i].setPosition(vecPos);
			rects[i].setSize(vecSize);
			*/
		}

		//
		// Do the things with the decibel meter
		//
		avgDbLevel = avg_history(&hisDbLevel);
		int h = (int)(win_wi * norm(avgDbLevel, 0, 90));
		int s = (int)lerp(prevAvgDbLevel, h, 0.8f);
		prevAvgDbLevel = h;
		beatRect.setSize(sf::Vector2f(s, 5));
		
		if (isSilence) 
		{
			beatRect.setSize(sf::Vector2f(win_wi, 5));
			beatRect.setFillColor(sf::Color::Red);
		} 
		else 
		{
			beatRect.setFillColor(sf::Color::Green);
		}

		//
		// Draw the beat bands
		//
		window.clear(backgroundColor); //sf::Color::Black);
		for (auto& r : rects) window.draw(r);

		window.draw(beatRect);
		window.display();

		Sleep(actualSleepTime); //buffersPerSec);
	}

	// cleanup
	printf("Cleaning up...\n");

	// probably forgot something here

	alcCaptureStop(pCaptureDevice);
	alcCloseDevice(pCaptureDevice);

	if ((err = alGetError()) != AL_NO_ERROR)
	{
		printf("Error code: %d", err);
	}

	fftwf_destroy_plan(plan);
	//fftwf_free(in);
	fftwf_free(out);

	in = NULL;
	out = NULL;
	plan = NULL;

	return 0;
}
