#pragma once

#define _USE_MATH_DEFINES

#include "fftw3.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <algorithm>


enum FFTSize : int
{
	FFT512 = 512,
	FFT1024 = 1024,
	FFT2048 = 2048,
	FFT4096 = 4096
};

class FFT
{
public:
	FFT(FFTSize size);
	~FFT();

	void process(float* outDataBuffer, float* in, size_t inLen);

	inline int getBufferSize() const { return size; }
	// 2 8-bit samples -> 1 float
	// which means the length is half the size
	inline int getFFTSize() const { return (size / 2); }
	inline int getMaxFftIndex() const { return (getFFTSize() / 2) + 1; }

	// in-place window functions
	static void winBlackman(float in[], size_t sampleLen);
	static void winHann(float in[], size_t sampleLen);
	static void winHamming(float in[], size_t sampleLen);
	static void winHarris(float in[], size_t sampleLen);

private:
	fftwf_plan plan;
	fftwf_complex* outData;
	float* intermediate;

	FFTSize size;
};
