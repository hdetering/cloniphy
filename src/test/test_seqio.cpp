#include <boost/test/unit_test.hpp>
#include <boost/timer/timer.hpp>

#include "../core/random.hpp"
#include "../core/seqio.hpp"
#include <map>
#include <memory>
#include <vector>
using namespace std;
using namespace seqio;

struct FixtureSeqio {
  FixtureSeqio() {
    BOOST_TEST_MESSAGE( "set up fixure" );
  }
  ~FixtureSeqio() {
    BOOST_TEST_MESSAGE( "teardown fixture" );
  }

  GenomeReference ref_genome;
};

BOOST_FIXTURE_TEST_SUITE( seqio, FixtureSeqio )

/* read reference genome from FASTA file */
BOOST_AUTO_TEST_CASE( read )
{
  string fn_fasta = "data/ref/Hg18.chr20-22.fa";

  BOOST_TEST_MESSAGE( "Reading genome from file: " << fn_fasta );
  ref_genome = GenomeReference(fn_fasta.c_str());
  BOOST_TEST_MESSAGE( "Read " << ref_genome.chromosomes.size() << " chromosomes." );
  BOOST_TEST_MESSAGE( "Indexing sequence records..." );
  ref_genome.indexRecords();

  BOOST_CHECK( ref_genome.chromosomes.size() == 3 );
}

/* generate reference genome */
BOOST_AUTO_TEST_CASE( gen )
{
  boost::timer::auto_cpu_timer t;
  //unsigned long len = 1000000; // genome length
  unsigned short num_chr = 10; // number of chromosomes
  unsigned long frag_len_mean = 100000; // mean chromosome length
  unsigned long frag_len_sd = 10000;
  vector<double> nuc_freqs = {0.3, 0.2, 0.2, 0.3};
  long seed = 123456789;
  RandomNumberGenerator<> rng(seed);

  BOOST_TEST_MESSAGE( "generating genome in " << num_chr << " fragments with length " << frag_len_mean << " (sd=" << frag_len_sd << ")..." );
  ref_genome = GenomeReference();
  ref_genome.generate(num_chr, frag_len_mean, frag_len_sd, nuc_freqs, rng);
  BOOST_TEST_MESSAGE( "done." );
  //BOOST_CHECK( ref_genome.length == len );
  BOOST_CHECK( ref_genome.num_records == num_chr );
  BOOST_CHECK( ref_genome.records.size() == num_chr );

  // write genome to file
  BOOST_TEST_MESSAGE( "writing generated genome to file 'ref_genome.fa'..." );
  ofstream f_fasta;
  f_fasta.open("ref.fa");
  writeFasta(ref_genome.records, f_fasta);
  f_fasta.close();
  BOOST_TEST_MESSAGE( "done." );

  // display summary stats
  ref_genome.indexRecords();
  BOOST_TEST_MESSAGE( "Genomic 3-mer counts:" );
  for (auto kv : ref_genome.map_3mer_pos) {
    BOOST_TEST_MESSAGE( "  " << kv.first << ": " << kv.second.size() );
  }
}

/* generate reference exome */
BOOST_AUTO_TEST_CASE( exome )
{
  boost::timer::auto_cpu_timer t;
  //unsigned long len = 1000000; // genome length
  unsigned int num_chr = 20000; // number of chromosomes
  unsigned long frag_len_mean = 300; // mean chromosome length
  unsigned long frag_len_sd = 200;
  vector<double> nuc_freqs = {0.25, 0.25, 0.25, 0.25};
  long seed = 123456789;
  RandomNumberGenerator<> rng(seed);

  BOOST_TEST_MESSAGE( "generating genome in " << num_chr << " fragments with length " << frag_len_mean << " (sd=" << frag_len_sd << ")..." );
  ref_genome = GenomeReference();
  ref_genome.generate(num_chr, frag_len_mean, frag_len_sd, nuc_freqs, rng);
  BOOST_TEST_MESSAGE( "done." );
  //BOOST_CHECK( ref_genome.length == len );
  BOOST_CHECK( ref_genome.num_records == num_chr );
  BOOST_CHECK( ref_genome.records.size() == num_chr );

  // write genome to file
  BOOST_TEST_MESSAGE( "writing generated exome to file 'exome.fa'..." );
  ofstream f_fasta;
  f_fasta.open("exome.fa");
  writeFasta(ref_genome.records, f_fasta);
  f_fasta.close();
  BOOST_TEST_MESSAGE( "done." );

  // display summary stats
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

BOOST_AUTO_TEST_CASE( tmap )
{
  string fn_fasta = "data/ref/min.fa";
  ref_genome = GenomeReference(fn_fasta.c_str());

  map<int, GenomeInstance> m;

  m[0] = GenomeInstance(ref_genome);
  m[1] = m[0];
  BOOST_TEST_MESSAGE( "m[0].vec_chr_len.size() : " << m[0].vec_chr_len.size() );
  BOOST_TEST_MESSAGE( "m[1].vec_chr_len.size() : " << m[1].vec_chr_len.size() );
  // BOOST_TEST_MESSAGE( "  m[1].vec_chr_len.push_back(42);" );
  // m[1].vec_chr_len.push_back(42);
  BOOST_TEST_MESSAGE( "  GenomeInstance gi = m[i];\n  gi.vec_chr_len.push_back(42);" );
  GenomeInstance gi = m[1];
  gi.vec_chr_len.push_back(42);
  m[1] = gi;
  BOOST_TEST_MESSAGE( "m[0].vec_chr_len.size() : " << m[0].vec_chr_len.size() );
  BOOST_TEST_MESSAGE( "m[1].vec_chr_len.size() : " << m[1].vec_chr_len.size() );
}

BOOST_AUTO_TEST_SUITE_END()
