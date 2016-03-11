#pragma once

#include "IAudioSource.h"
#include "FFT.h"

#include "al.h"
#include "alc.h"
#include <stdio.h>

class OpenALSource : IAudioSource
{
public:
	OpenALSource(int bufferSize, unsigned int freq);
	~OpenALSource();

	void getSample(unsigned char* buffer);
	
	void start();
	void stop();

	const ALCchar defaultDevice() const { return *pDeviceDefault; }

private:
	ALCint dataLen;
	ALCenum err;
	ALCuint freq;
	ALCsizei bufferSize;

	const ALCchar *pDeviceDefault;
	ALCdevice *pCaptureDevice;
};

