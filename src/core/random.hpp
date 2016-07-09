#ifndef RANDOM_H
#define RANDOM_H

#include <algorithm>
#include <functional>
#include <boost/bind.hpp>
#include <random>
#include <vector>
// Choosing the random number generator. (mt19937: Mersenne-Twister)
//typedef boost::mt19937 base_generator_type;
typedef std::mt19937 base_generator_type;

template <typename GeneratorType = base_generator_type>
struct RandomNumberGenerator {
	GeneratorType generator;

	//RandomNumberGenerator(long seed);
	RandomNumberGenerator(long seed) {
	  generator.seed(seed);
	}

	std::function<double()> getRandomFunctionDouble(double min, double max) {
		std::uniform_real_distribution<> dist(min, max);
		return boost::bind(dist, boost::ref(generator));
	}

	template <typename T>
	std::function<T()> getRandomFunctionInt(T min, T max) {
		std::uniform_int_distribution<> dist(min, max);
		return boost::bind(dist, boost::ref(generator));
	}

	template <typename IndexType = int>
	std::function<IndexType()> getRandomIndexWeighted(std::vector<double> weights) {
		std::discrete_distribution<IndexType> dist(weights.begin(), weights.end());
		return boost::bind(dist, boost::ref(generator));
	}

	std::function<double()> getRandomGamma(double shape, double scale) {
		std::gamma_distribution<> dist(shape, scale);
		return boost::bind(dist, boost::ref(generator));
	}

	template <typename RealType = double>
	std::function<RealType()> getRandomGammaMeanSd(RealType mean, RealType sd) {
		//std::uniform_real_distribution<> dist(min, max);
		//return boost::bind(dist, boost::ref(generator));
		// calculate gamma parameters
		double shape = float(mean*mean) / float(sd*sd);
		double scale = float(sd*sd) / float(mean);
		return getRandomGamma(shape, scale);
	}

  /** algorithm proposed in:
	 *  Rubin, Donald B. The Bayesian Bootstrap.
	 *  Ann. Statist. 9 (1981), no. 1, 130--134. doi:10.1214/aos/1176345338.
   */
	std::vector<double> getRandomProbs(int n) {
		std::vector<double> p(n);
		std::vector<double> r(n+1);
		r[0] = 0.0;
		r[1] = 1.0;
		auto random_double = getRandomFunctionDouble(0.0, 1.0);
		for (int i=2; i<n+1; ++i) {
			r[i] = random_double();
		}
		std::sort(r.begin(), r.end());
		for (int i=0; i<n; ++i) {
			p[i] = r[i+1] - r[i];
		}

		return p;
	}
};

/** Selects random element from container */
template <typename RandomGenerator = base_generator_type>
struct random_selector
{
	//On most platforms, you probably want to use std::random_device("/dev/urandom")()
	random_selector(RandomGenerator g = RandomGenerator(std::random_device("/dev/urandom")()))
		: _gen(g) {}

	// provide a seed to initialize the random generator
	random_selector(long seed) : _gen(RandomGenerator(seed)) {}

	template <typename Iter>
	Iter select(Iter start, Iter end) {
		std::uniform_int_distribution<> dist(0, std::distance(start, end) - 1);
		std::advance(start, dist(_gen));
		return start;
	}

	//convenience function
	template <typename Iter>
	Iter operator()(Iter start, Iter end) {
		return select(start, end);
	}

	//convenience function that works on anything with a sensible begin() and end(),
	// and returns with a ref to the value type
	template <typename Container>
	auto operator()(const Container& c) -> decltype(*begin(c))& {
		return *select(begin(c), end(c));
	}

private:
	RandomGenerator _gen;
};

#endif /* RANDOM_H */
