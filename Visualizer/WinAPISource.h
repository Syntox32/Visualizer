#pragma once

#include "IAudioSource.h"

class WinAPISource : public IAudioSource
{
public:
	WinAPISource();
	~WinAPISource();

	void getSample(unsigned char* buffer);
	const char* getDevice() const;
};

