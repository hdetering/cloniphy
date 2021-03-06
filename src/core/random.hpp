#ifndef RANDOM_H
#define RANDOM_H

#include "pcg-cpp/pcg_random.hpp" // pcg32, pcg_extras::seed_seq_from()

#include <algorithm> // std::sort()
#include <cassert>
#include <cmath> // pow()
#include <functional> // std::function<>, std::bind(), std::ref()
#include <random>
#include <vector>
// Choosing the random number generator. (mt19937: Mersenne-Twister)
//typedef std::mt19937 base_generator_type;
typedef pcg32 base_generator_type;
using std::pow;
// defined in <functional>
using std::bind;
using std::ref;

//template <typename GeneratorType = base_generator_type>
struct RandomNumberGenerator {

  // GeneratorType generator;
  base_generator_type generator;
  // Function generates numbers between zero and one (frequently needed).
  std::function<double()> rand_dbl = this->getRandomFunctionReal(0.0, 1.0);

  //RandomNumberGenerator(long seed);
  RandomNumberGenerator(long seed) {
    generator.seed(seed);
  }

  template <typename RealType = double>
  std::function<RealType()> 
  getRandomFunctionReal (
    const RealType min, 
    const RealType max
  ) 
  {
	  assert( max > min );
	  std::uniform_real_distribution<> dist(min, max);
	  return bind(dist, ref(generator));
  }

  template <typename IntType = int>
  std::function<IntType()> 
  getRandomFunctionInt(
    IntType min, 
    IntType max
  ) 
  {
	  std::uniform_int_distribution<> dist(min, max);
	  return bind(dist, ref(generator));
  }

  /** Returns a function that picks a random index using a set of weights. */
  template <typename IndexType = int, typename WeightType = double>
  std::function<IndexType()> 
  getRandomIndexWeighted (
    const std::vector<WeightType> weights
  ) 
  {
	  std::discrete_distribution<IndexType> dist(weights.begin(), weights.end());
	  return bind(dist, ref(generator));
  }

  /**
   * Get Poisson-distributed random function.
   * \param mean  Mean and variance. 
   */
  template <typename IntType = int, typename RealType = double>
  std::function<IntType()> 
  getRandomFunctionPoisson (
    const RealType mean
  ) 
  {
	  std::poisson_distribution<> dist(mean);
	  return bind(dist, ref(generator));
  }

  /**
   * Get Gamma-distributed random function. 
   */
  template <typename RealType = double>
  std::function<RealType()> 
  getRandomFunctionGamma (
    const RealType shape, 
    const RealType scale
  )
  {
	  std::gamma_distribution<> dist(shape, scale);
	  return bind(dist, ref(generator));
  }

  template <typename RealType = double>
  std::function<RealType()> 
  getRandomFunctionGammaMeanSd (
    const RealType mean, 
    const RealType sd
  )
  {
    //std::uniform_real_distribution<> dist(min, max);
    //return boost::bind(dist, boost::ref(generator));
    if (sd > 0) {
      // calculate gamma parameters
      double shape = float(mean*mean) / float(sd*sd);
      double scale = float(sd*sd) / float(mean);
      return getRandomFunctionGamma(shape, scale);
    } else {
      return [mean] () { return mean; };
    }
  }

  template <typename RealType = double>
  std::function<RealType()> 
  getRandomFunctionExponential (
    const RealType lambda
  ) 
  {
	  std::exponential_distribution<RealType> dist(lambda);
	  return bind(dist, ref(generator));
  }

  /**
   * Returns a function that samples from the binomial dristribution.
   * 
   * \param n  number of trials
   * \param p  probability of success
   * \returns  random function sampling from Binom(n, k) 
   */
  template <typename IntType = int, typename RealType = double>
  std::function<IntType()>
  getRandomFunctionBinomial (
    const IntType n,
    const RealType p
  )
  {
    std::binomial_distribution<> dist(n, p);
    return bind(dist, ref(generator));
  }

  /** algorithm proposed in:
   *  Rubin, Donald B. The Bayesian Bootstrap.
   *  Ann. Statist. 9 (1981), no. 1, 130--134. doi:10.1214/aos/1176345338.
   */
  template <typename RealType = double, typename IntType = int>
  std::vector<RealType> 
  getRandomProbs (
    const IntType n
  ) 
  {
    std::vector<RealType> p(n);
    std::vector<RealType> r(n+1);
    r[0] = 0.0;
    r[1] = 1.0;
    std::function<RealType()> r_unif = getRandomFunctionReal(0.0, 1.0);
    for (int i=2; i<n+1; ++i) {
      r[i] = r_unif();
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
  template <typename RealType = double>
  std::vector<RealType> 
  getRandomDirichlet (
	  const std::vector<RealType> alpha
  ) 
  {
    std::vector<RealType> res;
    RealType cumsum = 0.0;
    for (RealType a : alpha) {
      std::function<RealType()> rgamma = getRandomFunctionGamma(a, 1.0);
      RealType x = rgamma();
      res.push_back(x);
      cumsum += x;
    }
	  for (RealType& q : res)
	    q /= cumsum;

	  return res;
  }

  /*  
   *   Random variates from the negative binomial distribution.
   *
   *  NOTES
   *
   *    x = the number of failures before the n-th success
   *
   *  CREDIT
   * 
   *    Code adapted from an original version kindly provided by: David Posada.
   * 
   *  REFERENCE
   *
   *    Devroye, L. (1986).
   *    Non-Uniform Random Variate Generation.
   *    New York:Springer-Verlag.  Pages 488 and 543.
   *
   *  METHOD
   *
   *    Generate lambda as gamma with shape parameter "dispersion" (aka size) and scale
   *    parameter "mean/dispersion".  Return a Poisson deviate with mean lambda.

  **** NOTE: Extracted from rnbinom.c R code:
      rpois(rgamma(size, (1 - prob) / prob));
      rpois(rgamma(size, mu / size));

    The negative binomial distribution with dispersion = n and prob = p has density

      p(x) = Gamma(x+n)/(Gamma(n) x!) p^n (1-p)^x

    for x = 0, 1, 2, ..., n > 0 and 0 < p <= 1.

    A negative binomial distribution can arise as a mixture of Poisson distributions 
    with mean distributed as a Γ (pgamma) distribution with scale parameter 
    (1 - prob)/prob and shape parameter dispersion. In this model prob = scale/(1+scale), 
    and the mean is dispersion * (1 - prob)/prob. The variance in this parametrization is n (1-p)/p^2.

    The alternative parameterization, often used in ecology, and the one used here, 
    is by the mean mu, and the dispersion parameter, where prob = dispersion/(dispersion+mu). 
    The variance is mu + mu^2/dispersion in this parametrization.
  */

  template <typename IntType = int, typename RealType = double>
  IntType
  getRandomNegativeBinomial (
    const RealType mean, 
    const RealType dispersion)
  {
    // parameters for gamma distribution
    RealType shape  = dispersion;
    RealType scale  = (dispersion+mean)/dispersion-1;

    // sample NB as Poisson with Gamma-distributed mean
    RealType rgamma = getRandomFunctionGamma(shape, scale)();
    RealType rpois  = getRandomFunctionPoisson(rgamma)();

    return rpois;
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
  template <typename RealType = double>
  RealType
  getRandomParetoBounded (
    const RealType a, 
    const RealType l, 
    const RealType h
  ) 
  {
	  assert( a > 0 );
	  assert( l > 0 );
	  assert( h > l );
	  std::function<RealType()> r_unif = getRandomFunctionReal(0, 1);
	  RealType u = r_unif();

	  RealType x = pow(-(u*pow(h,a)-u*pow(l,a)-pow(h,a))/(pow(h,a)*pow(l,a)),-1.0/a);
	  return x;
  }

  /** 
   * Returns a random value sampled from a power-law distribution.
   * NOTE: Formula does not produce the expected results!
   *  
   *  Transformation from a uniform variable \in [0.0,1.0] according to
   *  Wolfram: http://mathworld.wolfram.com/RandomNumber.html
   *
   *  \param min minimum value
   *  \param max maximum value
   *  \param r rate parameter (sensu P(x)=1/L^r)
   */
  template <typename IntType = unsigned long, typename RealType = double>
  IntType
  getRandomPowerLaw (
	  const IntType min,
	  const IntType max,
	  const RealType r
  )
  {
	  std::function<RealType()> r_unif = getRandomFunctionReal(0.0, 1.0);
	  RealType y = r_unif();
	  IntType x = pow(pow(min,r+1) - pow(min,r+1)*y + pow(min,r+1), 1/(r+1));
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
