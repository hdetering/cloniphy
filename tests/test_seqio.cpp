#include <boost/test/unit_test.hpp>
#include <boost/timer/timer.hpp>

#include "seqio.hpp"
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

/* check if benchmark genome has been indexed correctly */
BOOST_AUTO_TEST_CASE( index_genome )
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
