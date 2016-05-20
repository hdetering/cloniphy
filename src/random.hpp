#ifndef RANDOM_H
#define RANDOM_H

#include <boost/bind.hpp>
#include <boost/function.hpp>
#include <boost/random.hpp>
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
	  boost::uniform_real<> uni_dist(min, max);
	  boost::variate_generator<base_generator_type&, boost::uniform_real<> > uni(generator, uni_dist);

	  boost::function<double()> f;
	  f = boost::bind(uni_dist, generator);
	  return f;
	}

	boost::function<int()> getRandomFunctionInt(int min, int max) {
	  boost::uniform_int<> uni_dist(min, max);
	  boost::variate_generator<base_generator_type&, boost::uniform_int<> > uni(generator, uni_dist);

	  boost::function<int()> f;
	  f = boost::bind(uni_dist, generator);
	  return f;
	}
};

/** Selects random element from container */
template <typename RandomGenerator = base_generator_type>
struct random_selector
{
	//On most platforms, you probably want to use std::random_device("/dev/urandom")()
	random_selector(RandomGenerator g = RandomGenerator(std::random_device("/dev/urandom")()))
		: gen(g) {}

	// provide a seed to initialize the random generator
	random_selector(long seed)	: gen(RandomGenerator(seed)) {}

	template <typename Iter>
	Iter select(Iter start, Iter end) {
		boost::uniform_int<> dist(0, std::distance(start, end) - 1);
		std::advance(start, dist(gen));
		return start;
	}

	//convenience function
	template <typename Iter>
	Iter operator()(Iter start, Iter end) {
		return select(start, end);
	}

	//convenience function that works on anything with a sensible begin() and end(), and returns with a ref to the value type
	template <typename Container>
	auto operator()(const Container& c) -> decltype(*begin(c))& {
		return *select(begin(c), end(c));
	}

private:
	RandomGenerator gen;
};

#endif /* RANDOM_H */
