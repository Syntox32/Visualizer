#include "Visualizer.h"


Visualizer::Visualizer(const FFTSize size)
	: useSerial(true), size(size), freq(44100)
{
	portName = "COM3"; /* \\\\.\\ */

	// use std::unique_ptr here?
	source = new OpenALSource(size, freq);
	serial = new Serial((char*)portName);
	fft = new FFT(size);
	
	sleepTime = 30L; // ms
	dbMin = 0;
	dbMax = -90;
	freqMin = 20;
	freqMax = 2000;

	if (!serial->IsConnected()) {
		printf("Serial disabled: could not connect to serial port.\n");
		useSerial = false;
	}
	else {
		printf("Serial enabled: connection opened successfully.\n");
	}

	printInfo();
}

Visualizer::~Visualizer()
{ }

float Visualizer::getDbLevel(float in[], size_t len)
{
	float lvl = 0.0f;
	for (unsigned int i = 0; i < len; i++)
		lvl += 20.0f * log10(in[i]);

	lvl /= (float)len;

	// clamp db level to range 0 to 90
	return max(min(lvl, dbMin), dbMax) - dbMax;
}

void Visualizer::scaleFft(float in[], size_t inLen)
{
	// perform scaling, either logarithmically or exponentially
	for (unsigned int i = 0; i < inLen; i++)
	{
		//in[i] *= 20;
		//in[i] = 20.0f * log10(in[i]);
		in[i] = sqrt(in[i]) * 2;
		//in[i] = sqrtf(in[i]);
		//in[i] = 10.0f * log10(in[i]);
	}
}

void Visualizer::init()
{
	LinearColumnSpectrum spectrum;

	spectrum.init(this);
	spectrum.start();
}

void Visualizer::printInfo() const
{
	printf("Using default recording device:\n\t%s\n\n", source->getDevice());
	printf("Sleep time: %d ms\n", sleepTime);
	printf("Operating at frequency: %d\n", freq);
	printf("Format: AL_FORMAT_MONO16\n");
	printf("Capture buffer size: %d\n", fft->getBufferSize());
	printf("FFT buffer size: %d\n", fft->getFFTSize());

	printf("\nNumber of bands: %d\n", 20);
	printf("Min db: %d\n", dbMin);
	printf("Max db: %d\n", dbMax);
	printf("Min freq: %d\n", freqMin);
	printf("Max freq: %d\n", freqMax);
}

std::vector<Freq> Visualizer::genLinFreqLimits(int n, int min, int max)
{
	const int delta = max - min;
	const int freqDelta = (int)(delta / n);
	std::vector<Freq> freqs;

	for (unsigned int i = 0; i < n; i++)
	{
		Freq f;

		f.min = floorf(min + (i * freqDelta));
		f.max = floorf(min + ((i + 1) * freqDelta));

		freqs.push_back(f);
	}

	return freqs;
}

std::vector<Freq> Visualizer::genLogFreqLimits(int n, size_t sampleRate)
{
	std::vector<Freq> freqs;

	for (unsigned int i = 0; i < n; i++)
	{
		float avg = 0.0f;
		int lowFreq = 0;

		if (i == 0)
			lowFreq = 0;
		else
			lowFreq = (int)((sampleRate / 2) / (float)pow(2, n - i));

		int hiFreq = (int)((sampleRate / 2) / (float)pow(2, (n - 1) - i));

		int lowBound = lowFreq; //freq_to_index(sample_rate, lowFreq, n);
		int highBound = hiFreq; //freq_to_index(sample_rate, hiFreq, n);

		Freq f;

		f.min = lowBound;
		f.max = highBound;

		freqs.push_back(f);
	}

	return freqs;
}

std::vector<Freq> Visualizer::genExpFreqLimits(int n, int min, int max, size_t sampleRate)
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