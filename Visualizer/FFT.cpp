#include "FFT.h"


FFT::FFT(FFTSize size)
	: size(size)
{
	// www.fftw.org/doc/One_002dDimensional-DFTs-of-Real-Data.html#One_002dDimensional-DFTs-of-Real-Data
	// in = (fftwf_complex*)fftw_alloc_complex(sizeof(fftwf_complex) * fftSize);
	inData = new float[size];
	outData = (fftwf_complex*)fftw_alloc_complex(sizeof(fftwf_complex) * size);
	plan = fftwf_plan_dft_r2c_1d(size, inData, outData, FFTW_ESTIMATE);
}

FFT::~FFT()
{
	fftwf_destroy_plan(plan);
	fftwf_free(outData);
	delete inData;
}

int FFT::feed(float in[], size_t sampleLen)
{
	for (unsigned int i = 0; i < sampleLen; i++) 
	{
		inData[i] = in[i];
	}
}

int FFT::process()
{
	fftwf_execute(plan);
}

void FFT::winBlackman(float in[], size_t sampleLen)
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

void FFT::winHann(float in[], size_t sampleLen)
{
	for (unsigned int i = 0; i < sampleLen; i++)
	{
		float mult = 0.5f * (1.0f - cosf(2.0f * M_PI * i / (sampleLen - 1)));
		in[i] *= mult;
	}
}

void FFT::winHamming(float in[], size_t sampleLen)
{
	float a = 0.54f;
	float b = 0.46f;

	for (unsigned int i = 0; i < sampleLen; i++)
	{
		float mult = a - (b * cosf(2.0f * M_PI * i / (sampleLen - 1)));
		in[i] *= mult;
	}
}

void FFT::winHarris(float in[], size_t sampleLen)
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