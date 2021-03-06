#pragma once

#include <memory>
#include <random>

#include "Visualizer.h"
#include "Utils.h"
#include "Biquad.h"

#include <SFML\Graphics.hpp>
#include <SFML\Window.hpp>

class Visualizer;
class Utils;
struct Freq;

class LinearColumnSpectrum
{
public:
	LinearColumnSpectrum();
	~LinearColumnSpectrum();

	void init(const Visualizer*);
	void start();
	void stop();

	inline float delta() const { return deltaTime; }

private:
	void nextFrame();
	void applyFilter(float *inData, size_t inLen);

	unsigned char *intermediate;
	float *readyOutBuffer;
	float *fftDataBuffer;

	size_t bufferSize;
	size_t fftSize;
	size_t maxFftIndex;

	float dbLvl;
	float prevRectHeight;
	float *prevBandEnergies;
	History* hisDbLevel;
	History* hisAverage;
	size_t numBands;
	size_t maxHistoryEntries;
	size_t maxBandEntries;
	bool silent;

	std::vector<Freq> fExp;
	std::vector<Freq> fLog;
	std::vector<Freq> fLin;

	Biquad *lp; // lowpass
	Biquad *bp; // bandpass
	Biquad *pf; // peakfilter

	std::vector<History*> bandHistories;

	History *hisEnergy;
	float averageEnergy;
	unsigned int beatCount;
	sf::Clock beatClock;

	sf::Vector2f currPos;
	sf::Vector2f currSize;

	Visualizer* vz;
	sf::RenderWindow* window;
	sf::Color backgroundColor;
	std::vector<sf::RectangleShape*> rects;
	sf::Vector2f rectPos;
	sf::Vector2f rectSize;
	unsigned int winHi;
	unsigned int winWi;
	float rectWi;
	float deltaTime;
	sf::Clock deltaClock;
	sf::Time frameDelta;
	sf::RectangleShape dbMeter;


};

