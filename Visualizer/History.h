#pragma once

#include <list>

class History
{
public:
	History(size_t maxLen)
		: maxLen(maxLen)
	{ 
		history = new std::list<float>();
	}

	~History()
	{
		delete history;
	}

	void push(float val)
	{
		if (history->size() == maxLen)
			history->pop_back();
		
		history->push_front(val);
	}

	float sum()
	{
		float ret = 0.0f;
		
		for (auto& val : *history)
			ret += val;
		
		return ret;
	}

	float avg()
	{
		float sumVal = sum();
		return (float)(sumVal / history->size());
	}

	inline std::list<float> *getHistory() const { return history; }
	inline size_t getMaxLen() const { return maxLen; }

private:
	size_t maxLen;
	std::list<float> *history;
};

