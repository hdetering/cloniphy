#include <boost/test/unit_test.hpp>
#include <boost/timer/timer.hpp>

#include "../core/random.hpp"
#include <boost/format.hpp>
using boost::format;
using boost::str;
#include <fstream>

using namespace std;

struct FixtureRandom {
  FixtureRandom() {
    BOOST_TEST_MESSAGE( "setup fixure" );
  }
  ~FixtureRandom() {
    BOOST_TEST_MESSAGE( "teardown fixture" );
  }

  long seed = 123456789;
  RandomNumberGenerator<> gen = RandomNumberGenerator<>(seed);
};

BOOST_FIXTURE_TEST_SUITE( rng, FixtureRandom )

/* generate random gamma-distributed values */
BOOST_AUTO_TEST_CASE( gamma )
{
  double shape = 2.0;
  double scale = 0.5;
  function<double()> random_gamma = gen.getRandomGamma(shape, scale);
  unsigned num_values = 100000;
  vector<double> vec_values(num_values);

  BOOST_TEST_MESSAGE( "Simulating gamma-ditributed values." );
  BOOST_TEST_MESSAGE( str(format("  shape: %f") % shape) );
  BOOST_TEST_MESSAGE( str(format("  scale: %f") % scale) );
  for (auto i=0; i<num_values; ++i) {
    vec_values[i] = random_gamma();
  }

  auto fn = "random_gamma_k2.0_t0.25.csv";
  BOOST_TEST_MESSAGE( "Writing results to output file" );
  BOOST_TEST_MESSAGE( str(format("  %s") % fn) );
  ofstream fs;
  fs.open(fn);
  for (double v : vec_values) {
    fs << format("%f") % v << endl;
  }
  fs.close();
}

BOOST_AUTO_TEST_SUITE_END()
