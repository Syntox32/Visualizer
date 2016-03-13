#include "Utils.h"

float Utils::norm(float x, float min, float max)
{
	return (x - min) / (max - min);
}

float Utils::lerp(float start, float end, float t)
{
	return (1.0f - t) * start + t * end;
}

int Utils::freqToIndex(size_t sampleRate, float freq, size_t maxFftIndex)
{
	float ny = sampleRate / 2.0f;
	return (int)roundf((freq / ny) * maxFftIndex);
}

void Utils::shortToFloat(float *out, unsigned char *in, size_t len)
{
	for (unsigned int i = 0; i < len; i++)
	{
		int n = (in[i] << 8 | in[i + 1]);
		out[i / 2] = (float)n / 32768.0f; // signed short max value
	}
}

bool Utils::isSilence(History *his, int minDb)
{
	bool ret = true;

	for (auto& val : *his->getHistory())
	{
		if ((int)round(val) > minDb)
			ret = false;
	}

	return ret;
}

