#include "FFT.h"


FFT::FFT(FFTSize size)
	: size(size)
{
	// www.fftw.org/doc/One_002dDimensional-DFTs-of-Real-Data.html#One_002dDimensional-DFTs-of-Real-Data
	size_t fftSize = getFFTSize();

	intermediate = new float[fftSize];
	outData = (fftwf_complex*)fftw_alloc_complex(sizeof(fftwf_complex) * fftSize);
	plan = fftwf_plan_dft_r2c_1d(fftSize, intermediate, outData, FFTW_ESTIMATE);
}

FFT::~FFT()
{
	fftwf_destroy_plan(plan);
	fftwf_free(outData);
	delete intermediate;
}

void FFT::process(float *outDataBuffer, float* in, size_t inLen)
{
	for (unsigned int k = 0; k < inLen; k++) 
	{
		intermediate[k] = in[k];
	}

	fftwf_execute(plan);

	float mag = 0.0f;
	unsigned int index = getMaxFftIndex();

	for (unsigned int i = 0; i < index; i++)
	{
		mag = sqrtf((outData[i][0] * outData[i][0]) + (outData[i][1] * outData[i][1]));
		outDataBuffer[i] = mag;
	}
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
