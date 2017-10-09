#include <boost/test/unit_test.hpp>
#include <boost/timer/timer.hpp>

#include "../core/random.hpp"
#include <boost/format.hpp>
#include <chrono>
#include <fstream>
#include <iostream>

using boost::format;
using boost::str;
using namespace std;
using boost::format;
using boost::str;


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

BOOST_FIXTURE_TEST_SUITE( rng , FixtureRandom )

/* pick indices according to weights */
BOOST_AUTO_TEST_CASE( indices )
{
  BOOST_TEST_MESSAGE( " Picking indices based on weights (0.1, 0.2, 0,3, 0.4)... " );

  vector<int> counts(4, 0);
  vector<double> probs = { 0.1, 0.2, 0.3, 0.4 };
  function<int()> rand_idx = gen.getRandomIndexWeighted(probs);

  for (int i=0; i<100000; ++i) {
    counts[rand_idx()]++;
  }

  for (int i=0; i<4; ++i) {
    BOOST_TEST_MESSAGE( format("%d: %d") % i % counts[i] );
  }
}

/* generate random gamma-distributed values */
BOOST_AUTO_TEST_CASE( gamma )
{
  double shape = 2.0;
  double scale = 0.5;
  function<double()> random_gamma = gen.getRandomGamma(shape, scale);
  unsigned num_values = 100000;
  vector<double> vec_values(num_values);

  BOOST_TEST_MESSAGE( "Simulating gamma-distributed values." );
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

/* generate random Dirichlet-distributed values */
BOOST_AUTO_TEST_CASE( dirichlet )
{
  vector<double> a, q;
  BOOST_TEST_MESSAGE( "Simulating Dirichlet-distributed values..." );

  a = {0.1, 0.1, 0.1};
  BOOST_TEST_MESSAGE( str(format("\n  alpha: [ %.4f, %.4f, %.4f ]") % a[0] % a[1] % a[2]) );
  for (int i=0; i<10; i++) {
    q = gen.getRandomDirichlet(a);
    BOOST_TEST_MESSAGE( str(format("  %d: [ %.4g, %.4g, %.4g ]") % i % q[0] % q[1] % q[2]) );
  }

  a = {1, 1, 1};
  BOOST_TEST_MESSAGE( str(format("\n  alpha: [ %.4f, %.4f, %.4f ]") % a[0] % a[1] % a[2]) );
  for (int i=0; i<10; i++) {
    q = gen.getRandomDirichlet(a);
    BOOST_TEST_MESSAGE( str(format("  %d: [ %.4g, %.4g, %.4g ]") % i % q[0] % q[1] % q[2]) );
  }

  a = {10, 10, 10};
  BOOST_TEST_MESSAGE( str(format("\n  alpha: [ %.4f, %.4f, %.4f ]") % a[0] % a[1] % a[2]) );
  for (int i=0; i<10; i++) {
    q = gen.getRandomDirichlet(a);
    BOOST_TEST_MESSAGE( str(format("  %d: [ %.4g, %.4g, %.4g ]") % i % q[0] % q[1] % q[2]) );
  }

  a = {100, 100, 100};
  BOOST_TEST_MESSAGE( str(format("\n  alpha: [ %.4f, %.4f, %.4f ]") % a[0] % a[1] % a[2]) );
  for (int i=0; i<10; i++) {
    q = gen.getRandomDirichlet(a);
    BOOST_TEST_MESSAGE( str(format("  %d: [ %.4g, %.4g, %.4g ]") % i % q[0] % q[1] % q[2]) );
  }
}

/* generate random numbers following a power-law distribution */
BOOST_AUTO_TEST_CASE( power )
{
  unsigned long min = 1000;
  unsigned long max = 1000000;
  double rate = 1.0;

  vector<unsigned long> vec_values;
  for (int i=0; i<100000; ++i) {
    vec_values.push_back(gen.getRandomPowerLaw(min, max, rate));
  }

  auto fn = "random_pow_r1.csv";
  BOOST_TEST_MESSAGE( "Writing results to output file" );
  BOOST_TEST_MESSAGE( str(format("  %s") % fn) );
  ofstream fs;
  fs.open(fn);
  for (double v : vec_values) {
    fs << format("%f") % v << endl;
  }
  fs.close();
}

/* generate random numbers following a Bounded Pareto distribution */
BOOST_AUTO_TEST_CASE( pareto )
{
  unsigned long min = 1000;
  unsigned long max = 1000000;
  double shape = 1.0;

  vector<double> vec_values;
  for (int i=0; i<100000; ++i) {
    vec_values.push_back(gen.getRandomParetoBounded(shape, double(min)/max, 1));
  }

  auto fn = "random_bp_r1.csv";
  BOOST_TEST_MESSAGE( "Writing results to output file" );
  BOOST_TEST_MESSAGE( str(format("  %s") % fn) );
  ofstream fs;
  fs.open(fn);
  for (double v : vec_values) {
    fs << format("%f") % v << endl;
  }
  fs.close();
}

BOOST_AUTO_TEST_CASE ( clock )
{
  cout << "system_clock" << endl;
  cout << chrono::system_clock::period::num << endl;
  cout << chrono::system_clock::period::den << endl;
  cout << "steady = " << boolalpha << chrono::system_clock::is_steady << endl << endl;
 
  cout << "high_resolution_clock" << endl;
  cout << chrono::high_resolution_clock::period::num << endl;
  cout << chrono::high_resolution_clock::period::den << endl;
  cout << "steady = " << boolalpha << chrono::high_resolution_clock::is_steady << endl << endl;

  cout << "steady_clock" << endl;
  cout << chrono::steady_clock::period::num << endl;
  cout << chrono::steady_clock::period::den << endl;
  cout << "steady = " << boolalpha << chrono::steady_clock::is_steady << endl << endl;
}

BOOST_AUTO_TEST_SUITE_END()
