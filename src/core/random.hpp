#ifndef RANDOM_H
#define RANDOM_H

#include <algorithm>
#include <boost/bind.hpp>
#include <boost/function.hpp>
#include <boost/random.hpp>
#include <vector>
// Choosing the random number generator. (mt19937: Mersenne-Twister)
typedef boost::mt19937 base_generator_type;

template <typename GeneratorType = base_generator_type>
struct RandomNumberGenerator {
	GeneratorType generator;

	//RandomNumberGenerator(long seed);
	RandomNumberGenerator(long seed) {
	  generator.seed(seed);
	}

	boost::function<double()> getRandomFunctionDouble(double min, double max) {
	  boost::uniform_real<> dist(min, max);
	  boost::variate_generator<base_generator_type&, boost::uniform_real<> > uni(generator, dist);
	  boost::function<double()> f = [dist, this] { return dist(generator); };
	  return f;
	}

	template <typename T>
	boost::function<T()> getRandomFunctionInt(T min, T max) {
	  boost::uniform_int<> dist(min, max);
	  boost::variate_generator<base_generator_type&, boost::uniform_int<> > uni(generator, dist);
	  boost::function<T()> f = [dist, this] { return dist(generator); };
	  return f;
	}

	boost::function<int()> getRandomIndexWeighted(std::vector<double> probabilities) {
		boost::random::discrete_distribution<> dist(probabilities.begin(), probabilities.end());
		boost::function<int()> f = [dist, this] { return dist(generator); };
		return f;
	}

	boost::function<double()> getRandomGamma(double shape, double scale) {
		boost::random::gamma_distribution<> dist(shape, scale);
		boost::random::variate_generator<base_generator_type&, boost::random::gamma_distribution<>> rand_gamma(generator, dist);
		//boost::function<double()> f = [dist, this] { return dist(generator); };
		boost::function<double()> f = boost::bind(dist, boost::ref(generator));
		return f;
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
		boost::uniform_int<> dist(0, std::distance(start, end) - 1);
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
