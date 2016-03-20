#pragma once

#include "IAudioSource.h"
#include "FFT.h"

#include "al.h"
#include "alc.h"
#include <stdio.h>

class OpenALSource : public IAudioSource
{
public:
	OpenALSource(const int bufferSize, unsigned int freq);
	~OpenALSource();

	void getSample(unsigned char* buffer);
	
	void start();
	void stop();

	const char* getDevice() const;

private:
	ALCint dataLen;
	ALCenum err;
	ALCuint freq;
	const ALCsizei bufferSize;
	unsigned char* intermediate;

	const ALCchar *pDeviceDefault;
	ALCdevice *pCaptureDevice;
	
	bool started;
};

