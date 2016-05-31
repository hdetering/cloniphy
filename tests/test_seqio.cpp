#include <boost/test/unit_test.hpp>
#include <boost/timer/timer.hpp>

#include <vector>
#include "random.hpp"
#include "seqio.hpp"
using namespace std;
using namespace seqio;

struct FixtureSeqio {
  FixtureSeqio() : ref_genome( "data/Hg18.chr20-22.fa" ) {
    BOOST_TEST_MESSAGE( "set up fixure" );
  }
  ~FixtureSeqio() {
    BOOST_TEST_MESSAGE( "teardown fixture" );
  }

  Genome ref_genome;
};

BOOST_FIXTURE_TEST_SUITE( seqio, FixtureSeqio )

/* generate */
BOOST_AUTO_TEST_CASE( generate )
{
  boost::timer::auto_cpu_timer t;
  unsigned long len = 1000000; // genome length
  unsigned short num_chr = 2; // number of chromosomes
  vector<double> nuc_freqs = {0.3, 0.2, 0.2, 0.3};
  long seed = 123456789;
  RandomNumberGenerator<> rng(seed);

  BOOST_TEST_MESSAGE( "generating genome of length " << len << "..." );
  ref_genome = Genome();
  ref_genome.generate(len, num_chr, nuc_freqs, rng);
  BOOST_TEST_MESSAGE( "done." );
  BOOST_CHECK( ref_genome.length == len );
  BOOST_CHECK( ref_genome.num_records == num_chr );
  BOOST_CHECK( ref_genome.records.size() == num_chr );

  ref_genome.indexRecords();
}

/* check if benchmark genome has been indexed correctly */
BOOST_AUTO_TEST_CASE( index )
{
  boost::timer::auto_cpu_timer t;
  ref_genome.indexRecords();
  BOOST_CHECK( ref_genome.length == 159071719 );
  BOOST_CHECK( ref_genome.num_records == 3 );
  BOOST_CHECK( ref_genome.records.size() == 3 );
  BOOST_CHECK( ref_genome.nuc_pos[0].size() == 35749166 );
  BOOST_CHECK( ref_genome.nuc_pos[1].size() == 28451472 );
  BOOST_CHECK( ref_genome.nuc_pos[2].size() == 28496320 );
  BOOST_CHECK( ref_genome.nuc_pos[3].size() == 35829712 );
}

BOOST_AUTO_TEST_SUITE_END()
