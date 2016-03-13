#pragma once

#include "FFT.h"
#include "Utils.h"
#include "SerialClass.h"
#include "IAudioSource.h"
#include "OpenALSource.h"
#include "WinAPISource.h"
#include "SpectrumDisplay.h"

#include <vector>

struct Freq
{
	float min;
	float max;
};

class Visualizer
{
public:
	Visualizer(FFTSize size);
	~Visualizer();

	void init();

	inline unsigned int bufferSize() const { return (unsigned int)size; }
	inline unsigned int fftSize() const { return (unsigned int)(size / 2); }

	std::vector<Freq> genLinFreqLimits(int n, int min, int max);
	std::vector<Freq> genLogFreqLimits(int n, size_t sampleRate);
	std::vector<Freq> genExpFreqLimits(int n, int min, int max, size_t sampleRate);

private:
	SpectrumDisplay* display = nullptr;
	IAudioSource* source = nullptr;
	Serial* serial = nullptr;
	FFT* fft = nullptr;

	const char* portName;
	const unsigned int freq;
	const FFTSize size;
	bool useSerial;

	long sleepTime;
	int dbMin;
	int dbMax;
	int freqMin;
	int freqMax;

	void printInfo();
};

