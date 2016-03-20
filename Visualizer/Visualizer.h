#pragma once

#include "FFT.h"
#include "SerialClass.h"
#include "IAudioSource.h"
#include "OpenALSource.h"
#include "WinAPISource.h"
#include "LinearColumnSpectrum.h"

#include <vector>
#include <memory>

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

	inline unsigned int getBufferSize() const { return (unsigned int)size; }
	inline unsigned int getFFTSize() const { return (unsigned int)(size / 2); }

	std::vector<Freq> genLinFreqLimits(int n, int min, int max);
	std::vector<Freq> genLogFreqLimits(int n, size_t sampleRate);
	std::vector<Freq> genExpFreqLimits(int n, int min, int max, size_t sampleRate);

	void printInfo() const;
	float getDbLevel(float in[], size_t len); //, int minDb, int maxDb);
	void scaleFft(float in[], size_t inLen);

	IAudioSource* source;
	Serial* serial;
	FFT* fft;
	long sleepTime;
	bool useSerial;
	const FFTSize size;
	
	const unsigned int freq;
	const char* portName;

	int dbMin;
	int dbMax;
	int freqMin;
	int freqMax;
};

