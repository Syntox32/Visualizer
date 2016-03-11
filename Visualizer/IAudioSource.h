#pragma once

class IAudioSource
{
public:
	virtual unsigned char getSample(size_t width) = 0;
};
