#pragma once

#include <list>
#include "History.h"

class Utils
{
public:
	static float norm(float x, float min, float max);
	static float lerp(float start, float end, float t);

	static int freqToIndex(size_t sampleRate, float freq, size_t maxFftIndex);
	static void shortToFloat(float *out, unsigned char *in, size_t len);
	static bool isSilence(History *his, int minDb);

	//float average(float in[], size_t len);
	//float variance(float in[], size_t len, float energy);
};

