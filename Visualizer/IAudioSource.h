#pragma once

class IAudioSource
{
public:
	virtual void start() = 0;
	virtual void stop() = 0;

	virtual void getSample(unsigned char* buffer) = 0;
	virtual const char* getDevice() const = 0;
};
