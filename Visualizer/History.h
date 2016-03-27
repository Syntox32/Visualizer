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

	float var(float pivotValue)
	{
		float ret = 0.0f;

		for (auto& val : *history)
			ret += val - pivotValue;

		return (ret / (float)history->size());
	}

	// square average thing
	// something something statistics en.wikipedia.org/wiki/Mean_squared_error
	float sqAvg()
	{
		float ret = 0.0f;

		for (auto& val : *history)
			ret += val; //(val * val);

		return ((ret * ret) / (float)history->size());
	}

	// square variance
	float sqVar(float pivotValue)
	{
		float ret = 0.0f;

		for (auto& val : *history)
			ret += val - pivotValue; //powf(val - pivotValue, 2.0f);

		return ((ret * ret) / (float)history->size());
	}

	inline size_t count() const { return history->size(); }
	inline std::list<float> *getHistory() const { return history; }
	inline size_t getMaxLen() const { return maxLen; }

private:
	size_t maxLen;
	std::list<float> *history;
};

