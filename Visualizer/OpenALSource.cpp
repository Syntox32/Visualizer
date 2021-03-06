#include "OpenALSource.h"


OpenALSource::OpenALSource(const int bufferSize, unsigned int freq)
	: dataLen(0), freq(freq), bufferSize(bufferSize), started(false)
{
	intermediate = new unsigned char[bufferSize];

	pDeviceDefault = alcGetString(
		nullptr,
		ALC_CAPTURE_DEVICE_SPECIFIER); // get default capture device

	if ((err = alGetError()) != AL_NO_ERROR)
	{
		printf("Error code: %d\n", err);
	}

	// blitzbasic.com/Community/posts.php?topic=90830
	pCaptureDevice = alcCaptureOpenDevice(
		pDeviceDefault,
		freq,              // samples per second
		AL_FORMAT_MONO16,  // one channel with two 8-bit samples
		bufferSize);       // sample_rate * res * num_tracks * sec

	if ((err = alGetError()) != AL_NO_ERROR)
	{
		printf("Error code: %d\n", err);
	}

	if (!pCaptureDevice)
	{
		printf("Could not create capture device.\n");
	}
}

void OpenALSource::start() 
{
	printf("Starting capture device...\n");
	alcCaptureStart(pCaptureDevice);

	if ((err = alGetError()) != AL_NO_ERROR)
	{
		printf("Error code: %d\n", err);
	}
	
	started = true;
}

void OpenALSource::stop()
{
	alcCaptureStop(pCaptureDevice);

	if ((err = alGetError()) != AL_NO_ERROR)
	{
		printf("Error code: %d\n", err);
	}
}

const char* OpenALSource::getDevice() const
{
	return (const char*)pDeviceDefault;
}

void OpenALSource::getSample(unsigned char* buffer)
{
	if (!started)
		printf("Warning: Source not started.");

	// Get the samples
	alcGetIntegerv(
		pCaptureDevice,
		ALC_CAPTURE_SAMPLES,
		(ALCsizei)sizeof(ALCint),
		&dataLen);

	if ((err = alGetError()) != AL_NO_ERROR)
	{
		printf("Error code: %d\n", err);
	}


	// if we don't have enough samples yet,
	// we just wait another frame and just 
	// proceed with the previous one
	if (dataLen >= bufferSize)
	{
		alcCaptureSamples(
			pCaptureDevice,
			(ALCvoid*)buffer,
			bufferSize);

		//printf("dataLen: %d\n", dataLen);

		if ((err = alGetError()) != AL_NO_ERROR)
		{
			printf("Error code: %d\n", err);
		}
	}
	else
	{
		// we don't override the previous data in the buffer
		// and just use the previous values to to the next calculations

		//printf("Warning: dataLen < captureBufferLen %d\n", dataLen);
	}
}

OpenALSource::~OpenALSource()
{
	alcCloseDevice(pCaptureDevice);

	if ((err = alGetError()) != AL_NO_ERROR)
	{
		printf("Error code: %d\n", err);
	}
}
