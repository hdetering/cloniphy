#ifndef RANDOM_H
#define RANDOM_H

#include <algorithm>
#include <boost/bind.hpp>
#include <cassert>
#include <cmath> // pow()
#include <functional>
#include <random>
#include <vector>
// Choosing the random number generator. (mt19937: Mersenne-Twister)
//typedef boost::mt19937 base_generator_type;
typedef std::mt19937 base_generator_type;
using std::pow;

template <typename GeneratorType = base_generator_type>
struct RandomNumberGenerator {
	GeneratorType generator;

	//RandomNumberGenerator(long seed);
	RandomNumberGenerator(long seed) {
	  generator.seed(seed);
	}

	std::function<double()> getRandomFunctionDouble(double min, double max) {
		assert( max > min );
		std::uniform_real_distribution<> dist(min, max);
		return boost::bind(dist, boost::ref(generator));
	}

	template <typename T>
	std::function<T()> getRandomFunctionInt(T min, T max) {
		std::uniform_int_distribution<> dist(min, max);
		return boost::bind(dist, boost::ref(generator));
	}

	template <typename IndexType = int, typename WeightType = double>
	std::function<IndexType()> getRandomIndexWeighted(std::vector<WeightType> weights) {
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
		if (sd > 0) {
			// calculate gamma parameters
			double shape = float(mean*mean) / float(sd*sd);
			double scale = float(sd*sd) / float(mean);
			return getRandomGamma(shape, scale);
		} else {
			return [mean] () { return mean; };
		}

	}

	template <typename RealType = double>
	std::function<RealType()> getRandomExponential(double lambda) {
		std::exponential_distribution<RealType> dist(lambda);
		return boost::bind(dist, boost::ref(generator));
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

	/** algorithm proposed in:
	 *  Introduction to the Dirichlet distribution and related processes
   *  BA Frigyik, A Kapila, MR Gupta
   *  Dept. Elect. Eng., Univ. Washington, Seattle, WA, USA, UWEETR-2010-0006, 2010
	 */
	std::vector<double> getRandomDirichlet(std::vector<double> alpha) {
		std::vector<double> res;
		double cumsum = 0.0;
		for (double a : alpha) {
			auto rgamma = getRandomGamma(a, 1);
			double x = rgamma();
			res.push_back(x);
			cumsum += x;
		}
		for (double& q : res)
			q /= cumsum;

		return res;
	}

	/**
	 *  Internally samples from a uniform distribution (U~Unif(0,1)) and transforms to
	 *  Bounded Pareto by inverse-transform menthod:
	 *   x = \left(-\frac{U H^\alpha - U L^\alpha - H^\alpha}{H^\alpha L^\alpha}\right)^{-\frac{1}{\alpha}}
	 *
	 *  \param a shape
	 *  \param l minimum value (>0)
	 *  \param h maximum value (>l)
	 */
	double getRandomParetoBounded(const double a, const double l, const double h) {
		assert( a > 0 );
		assert( l > 0 );
		assert( h > l );
		std::function<double()> r_unif = getRandomFunctionDouble(0, 1);
		double u = r_unif();

		double x = pow(-(u*pow(h,a)-u*pow(l,a)-pow(h,a))/(pow(h,a)*pow(l,a)),-1.0/a);
		return x;
	}

	/** NOTE: Formula does not produce the expected results!
	 * Returns a random value sampled from a power-law distribution.
	 *  Transformation from a uniform variable \in [0.0,1.0] according to
	 *  Wolfram: http://mathworld.wolfram.com/RandomNumber.html
	 *
	 *  \param min minimum value
	 *  \param max maximum value
	 *  \param r rate parameter (sensu P(x)=1/L^r)
	 */
	unsigned long getRandomPowerLaw (
		const unsigned long min,
		const unsigned long max,
		const double r)
	{
		auto r_dbl_unif = getRandomFunctionDouble(0.0, 1.0);
		double y = r_dbl_unif();
		unsigned long x = pow(pow(min,r+1) - pow(min,r+1)*y + pow(min,r+1), 1/(r+1));
		return x;
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
		assert( c.size() > 0 );
		return *select(begin(c), end(c));
	}

private:
	RandomGenerator _gen;
};

#endif /* RANDOM_H */
