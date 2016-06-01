#include <boost/test/unit_test.hpp>

#include <boost/format.hpp>
#include <vector>
#include "yaml-cpp/yaml.h"

using namespace std;
using boost::format;

// forward declarations
template <typename T>
T get_param(const char* name, YAML::Node& config);

struct FixtureInterface {
  FixtureInterface() {
    BOOST_TEST_MESSAGE( "setup fixure" );
    config = YAML::LoadFile("config.yml");
  }
  ~FixtureInterface() {
    BOOST_TEST_MESSAGE( "teardown fixure" );
  }
  YAML::Node config = YAML::Node();
};

BOOST_FIXTURE_TEST_SUITE( cli, FixtureInterface )

BOOST_AUTO_TEST_CASE( params )
{
  BOOST_TEST_MESSAGE( "checking parameters in config file...\n" );

  BOOST_TEST_MESSAGE( "clones" << ": " << get_param<int>("clones", config) );
  BOOST_TEST_MESSAGE( "mutations" << ": " << get_param<int>("mutations", config) );
  BOOST_TEST_MESSAGE( "init-muts" << ": " << get_param<int>("init-muts", config) );
  vector<float> freqs = get_param<vector<float> >("freqs", config);
  BOOST_TEST_MESSAGE( format("freqs: [ %.2f, %.2f, %.2f, %.2f ]") % freqs[0] % freqs[1] % freqs[2] % freqs[3] );
  BOOST_TEST_MESSAGE( "ref-len" << ": " << get_param<unsigned long>("ref-len", config) );
  vector<double> rn_freqs = get_param<vector<double> >("ref-nuc-freqs", config);
  BOOST_TEST_MESSAGE( format("freqs: [ %.2f, %.2f, %.2f, %.2f ]") % rn_freqs[0] % rn_freqs[1] % rn_freqs[2] % rn_freqs[3] );
  BOOST_TEST_MESSAGE( "reference" << ": " << get_param<string>("reference", config) );
  BOOST_TEST_MESSAGE( "reference-vcf" << ": " << get_param<string>("reference-vcf", config) );
  BOOST_TEST_MESSAGE( "tree" << ": " << get_param<string>("tree", config) );
  BOOST_TEST_MESSAGE( "verbosity" << ": " << get_param<int>("verbosity", config) );
  BOOST_TEST_MESSAGE( "seed" << ": " << get_param<long>("seed", config) );
}

BOOST_AUTO_TEST_SUITE_END()

template <typename T>
T get_param(const char* name, YAML::Node& config) {
  T value;
  if (config[name])
    value = config[name].as<T>();
  return value;
}
