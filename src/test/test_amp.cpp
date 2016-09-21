#include <boost/test/unit_test.hpp>
#include <boost/timer/timer.hpp>

#include <boost/format.hpp>
using boost::format;
using boost::str;
#include <fstream>
#include "../core/amp.hpp"
#include "../core/seqio.hpp"
#include "../core/vario.hpp"
using namespace amp;

using namespace std;

struct FixtureAmp {
  FixtureAmp() {
    BOOST_TEST_MESSAGE( "setup fixure" );
  }
  ~FixtureAmp() {
    BOOST_TEST_MESSAGE( "teardown fixture" );
  }

  long seed = 123456789;
  RandomNumberGenerator<> rng = RandomNumberGenerator<>(seed);
};

BOOST_FIXTURE_TEST_SUITE( amp, FixtureAmp )

BOOST_AUTO_TEST_CASE( mda ) {
  boost::timer::auto_cpu_timer t;
  unsigned seq_len = 100000000;
  double amp_size_mean = 70000.0;
  double amp_size_sd = 30000.0;
  int fold = 1000;
  vector<unsigned> cvg;
  ofstream fs_cvg;

  BOOST_TEST_MESSAGE( "Simulating MDA..." );
  BOOST_TEST_MESSAGE( str(format("  template len:\t%d") % seq_len) );
  BOOST_TEST_MESSAGE( str(format("  amplicon len:\t%.0f (+/- %.0f)") % amp_size_mean % amp_size_sd) );
  BOOST_TEST_MESSAGE( str(format("  fold:\t%d") % fold) );

  BOOST_TEST_MESSAGE( "\nchromosome copy 1..." );
  simulateMda(cvg, seq_len, amp_size_mean, amp_size_sd, fold, rng);
  BOOST_TEST_MESSAGE( "  done." );
  BOOST_TEST_MESSAGE( "  writing results to file 'mda_cvg.dist.csv'..." );
  fs_cvg.open("mda_cvg.dist.csv");
  for (auto c : cvg)
    fs_cvg << c << endl;
  fs_cvg.close();
/*
  BOOST_TEST_MESSAGE( "\nchromosome copy 2." );
  simulateMda(cvg, seq_len, amp_size_mean, amp_size_sd, fold, rng);
  BOOST_TEST_MESSAGE( "  done." );
  BOOST_TEST_MESSAGE( "  writing results to file 'mda_cvg_2.csv'..." );
  fs_cvg.open("mda_cvg_2.csv");
  for (auto c : cvg)
    fs_cvg << c << endl;
  fs_cvg.close();
  */
}

BOOST_AUTO_TEST_CASE( art ) {
  boost::timer::auto_cpu_timer t;
  // reference genome
  int num_ref_seq = 1;
  unsigned seq_len = 100000000;
  unsigned seq_len_sd = 0;
  vector<double> nuc_freqs = { 0.3, 0.2, 0.2, 0.3 };
  // MDA params
  double amp_size_mean = 70000.0;
  double amp_size_sd = 30000.0;
  int fold = 1000;
  // local variables
  vector<unsigned> cvg;
  ofstream fs_cvg;

  BOOST_TEST_MESSAGE( "generating random genome (" << seq_len << " bp)... ");
  seqio::Genome ref_genome;
  ref_genome.generate(num_ref_seq, seq_len, seq_len_sd, nuc_freqs, rng);
  BOOST_TEST_MESSAGE( "  duplicating reference genome..." );
  ref_genome.duplicate();
  auto fn_ref = "ref.100M.fa";
  BOOST_TEST_MESSAGE( "  writing reference genome to file '" << fn_ref << "'...\n" );
  ofstream fs_ref;
  fs_ref.open(fn_ref);
  writeFasta(ref_genome.records, fs_ref);
  fs_ref.close();
  //BOOST_TEST_MESSAGE( "indexing genome... ");
  //ref_genome.indexRecords();

  BOOST_TEST_MESSAGE( "Simulating MDA..." );
  BOOST_TEST_MESSAGE( str(format("  template len:\t%d") % seq_len) );
  BOOST_TEST_MESSAGE( str(format("  amplicon len:\t%.0f (+/- %.0f)") % amp_size_mean % amp_size_sd) );
  BOOST_TEST_MESSAGE( str(format("  fold:\t%d") % fold) );

  string fn_cvg = "mda_cvg.dist.cfa";
  fs_cvg.open(fn_cvg);
  for (auto rec : ref_genome.records) {
    BOOST_TEST_MESSAGE( "\nchromosome '" << rec.id << "'..." );
    simulateMda(cvg, seq_len, amp_size_mean, amp_size_sd, fold, rng);
    //simulateMdaProcess(cvg, seq_len, amp_size_mean, amp_size_sd, fold, rng);
    BOOST_TEST_MESSAGE( "  done." );
    fs_cvg << ">" << rec.id << endl;
    for (auto c : cvg)
      fs_cvg << c << "#";
    fs_cvg << endl;
  }
  fs_cvg.close();
  BOOST_TEST_MESSAGE( "  wrote results to file '" << fn_cvg << "'" );
}

BOOST_AUTO_TEST_SUITE_END()
