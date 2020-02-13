#include <boost/test/unit_test.hpp>
#include <boost/timer/timer.hpp>

#include "../core/random.hpp"
#include "../core/seqio.hpp"
#include "../core/seqio/GenomeReference.hpp"
#include "../core/seqio/GenomeInstance.hpp"
#include <boost/icl/interval.hpp>
#include <boost/icl/interval_map.hpp>
#include <map>
#include <memory>
#include <set>
#include <vector>
using namespace std;
using namespace seqio;
using boost::icl::interval_map;
using boost::icl::interval;

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
  RandomNumberGenerator rng(seed);

  BOOST_TEST_MESSAGE( "generating genome in " << num_chr << " fragments with length " << frag_len_mean << " (sd=" << frag_len_sd << ")..." );
  ref_genome = GenomeReference();
  ref_genome.generate_nucfreqs(num_chr, frag_len_mean, frag_len_sd, nuc_freqs, rng);
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
  RandomNumberGenerator rng(seed);

  BOOST_TEST_MESSAGE( "generating genome in " << num_chr << " fragments with length " << frag_len_mean << " (sd=" << frag_len_sd << ")..." );
  ref_genome = GenomeReference();
  ref_genome.generate_nucfreqs(num_chr, frag_len_mean, frag_len_sd, nuc_freqs, rng);
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
BOOST_AUTO_TEST_CASE ( index )
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

BOOST_AUTO_TEST_CASE ( tmap )
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

/* tests related to the BOOST Interval Container Library (ICL)
 * http://www.boost.org/doc/libs/1_64_0/libs/icl/doc/html/index.html */
BOOST_AUTO_TEST_CASE ( icl ) {
  BOOST_TEST_MESSAGE ( "----------------------------------------------" );
  BOOST_TEST_MESSAGE ( "Testing Boost Interval Container Library (ICL)" );
  BOOST_TEST_MESSAGE ( "----------------------------------------------" );

  interval_map<unsigned long, int> imap_copies;

  // insert a couple of intervals

  auto i1 = interval<unsigned long>::right_open(1, 100);
  imap_copies += make_pair(i1, 1);
  BOOST_TEST_MESSAGE ( "\nadd: " << i1);

  auto i2 = interval<unsigned long>::right_open(80, 120);
  BOOST_TEST_MESSAGE ( "add: " << i2);
  imap_copies += make_pair(i2, 1);

  auto i3 = interval<unsigned long>::right_open(30, 60);
  BOOST_TEST_MESSAGE ( "add: " << i3);
  imap_copies += make_pair(i3, 1);

  BOOST_TEST_MESSAGE ( "\nresult:\n" << imap_copies << "\n" );

  BOOST_TEST_MESSAGE ( "\nintervals:" );
  for (auto const & i : imap_copies) {
    BOOST_TEST_MESSAGE ( "\t" << i.first << ", " << i.second );
  }
}

/* generate reference genome, introduce CNVs and write CN state-tiled output */
BOOST_AUTO_TEST_CASE( tile )
{
  boost::timer::auto_cpu_timer t;
  unsigned short num_chr = 3; // number of chromosomes
  unsigned long frag_len_mean = 100000; // mean chromosome length
  unsigned long frag_len_sd = 0;
  vector<double> nuc_freqs = {0.3, 0.2, 0.2, 0.3};
  long seed = 42;
  RandomNumberGenerator rng(seed);
  unsigned long start_bp, end_bp, len_bp = 0;

  BOOST_TEST_MESSAGE( "--------------------------------" );
  BOOST_TEST_MESSAGE( " Generate genome:" );
  BOOST_TEST_MESSAGE( "--------------------------------\n" );
  BOOST_TEST_MESSAGE( "chromosomes:\t" << num_chr );
  BOOST_TEST_MESSAGE( "chr length:\t" << frag_len_mean << " (+/- " << frag_len_sd << ")\n");
  ref_genome = GenomeReference();
  ref_genome.generate_nucfreqs(num_chr, frag_len_mean, frag_len_sd, nuc_freqs, rng);
  GenomeInstance genome(ref_genome);
  BOOST_TEST_MESSAGE( genome );

  BOOST_TEST_MESSAGE( "--------------------------------" );
  BOOST_TEST_MESSAGE( " Introduce CNVs:" );
  BOOST_TEST_MESSAGE( "--------------------------------\n" );

  string id_chr;
  shared_ptr<ChromosomeInstance> chr_copy;
  double start_rel = 0.0;
  double len_rel = 0.0;
  bool is_forward = true;
  bool is_telomeric = false;
  bool is_first_arm = false;

  BOOST_TEST_MESSAGE( "Delete second copy of chromosome 2" );
  BOOST_TEST_MESSAGE( "--------------------------------\n" );
  id_chr = "chr2";
  chr_copy = genome.map_id_chr[id_chr][1];
  genome.deleteChromosomeInstance( chr_copy );
  BOOST_TEST_MESSAGE( "result:\n" << genome );

  BOOST_TEST_MESSAGE( "" );
  BOOST_TEST_MESSAGE( "Amplify telomeric region on chromosome 1" );
  BOOST_TEST_MESSAGE( "-----------------------------------------\n" );
  id_chr = "chr1";
  chr_copy = genome.map_id_chr[id_chr][0];
  start_rel = 0.5;
  len_rel = 0.5;
  is_forward = true;
  is_telomeric = true;
  is_first_arm = false;
  chr_copy->amplifyRegion(
    start_bp, end_bp, len_bp,
    start_rel, 
    len_rel, 
    is_forward, 
    is_telomeric,
    is_first_arm
  );
  BOOST_TEST_MESSAGE( "result:\n" << genome );

  BOOST_TEST_MESSAGE( "Delete focal region on chromosome 3" );
  BOOST_TEST_MESSAGE( "-----------------------------------\n" );
  id_chr = "chr3";
  chr_copy = genome.map_id_chr[id_chr][1];
  start_rel = 0.25;
  len_rel = 0.5;
  is_forward = true;
  is_telomeric = false;
  is_first_arm = false;
  chr_copy->deleteRegion(
    start_bp, end_bp, len_bp,
    start_rel, 
    len_rel, 
    is_forward, 
    is_telomeric,
    is_first_arm
  );
  BOOST_TEST_MESSAGE( "result:\n" << genome );

  BOOST_TEST_MESSAGE( "Amplify focal region on chromosome 2" );
  BOOST_TEST_MESSAGE( "------------------------------------\n" );
  id_chr = "chr2";
  chr_copy = genome.map_id_chr[id_chr][0];
  start_rel = 0.75;
  len_rel = 0.5;
  is_forward = false;
  is_telomeric = false;
  is_first_arm = false;
  chr_copy->amplifyRegion(
    start_bp, end_bp, len_bp,
    start_rel, 
    len_rel, 
    is_forward, 
    is_telomeric,
    is_first_arm
  );
  BOOST_TEST_MESSAGE( "result:\n" << genome );

  BOOST_TEST_MESSAGE( "Indexing SegmentCopies" );
  BOOST_TEST_MESSAGE( "----------------------\n" );
  genome.indexSegmentCopies();
  BOOST_TEST_MESSAGE( "result:" );
  for (auto const & kv : genome.map_chr_seg) {
    string id_chr = kv.first;
    TSegMap imap_seg = kv.second;
    cout << "\t" << id_chr << endl;
    cout << "\t\t" << imap_seg << endl;
  }

  // NOTE: moved this functionality to bamio::BulkSampleGenerator
  // BOOST_TEST_MESSAGE( "--------------------------------" );
  // BOOST_TEST_MESSAGE( " Write tiled genome to files" );
  // BOOST_TEST_MESSAGE( "--------------------------------\n" );

  // string fn_pfx = "test.segio_tile";
  // genome.writeFastaTiled(ref_genome, fn_pfx, 10, 500);
}

BOOST_AUTO_TEST_SUITE_END()
