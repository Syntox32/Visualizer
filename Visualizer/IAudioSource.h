#pragma once

class IAudioSource
{
public:
	virtual void getSample(unsigned char* buffer) = 0;
	virtual const char* getDevice() const = 0;
};
