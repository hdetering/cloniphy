#include <boost/test/unit_test.hpp>
#include <boost/timer/timer.hpp>

#include "../core/random.hpp"
#include "../core/seqio.hpp"
#include "../core/vario.hpp"
#include <boost/format.hpp>
#include <boost/icl/interval_map.hpp>
using namespace boost::icl;
#include <fstream>
#include <vector>
using boost::format;
using evolution::SubstitutionModel;
using namespace std;
using namespace vario;

struct FixtureVario {
  FixtureVario() {
    BOOST_TEST_MESSAGE( "setup fixure" );
    fn_vcf = "data/germline.vcf";

    unsigned long ref_genome_len = 1000000;
    int num_ref_seq = 10;
    int seq_len_mean = 100000;
    int seq_len_sd = 10000;
    vector<double> nuc_freqs = { 0.3, 0.2, 0.2, 0.3 };
    // initialize random number generator

    BOOST_TEST_MESSAGE( "generating random genome (" << ref_genome_len << " bp)... ");
    ref_genome.generate(num_ref_seq, seq_len_mean, seq_len_sd, nuc_freqs, rng);
    auto fn_ref = "ref.fa";
    BOOST_TEST_MESSAGE( "writing reference genome to file '" << fn_ref << "'...\n" );
    ofstream fs_ref;
    fs_ref.open(fn_ref);
    writeFasta(ref_genome.records, fs_ref);
    fs_ref.close();

    BOOST_TEST_MESSAGE( "indexing genome... ");
    ref_genome.indexRecords();

    double Q[4][4] = {
      {    0.0, 0.0406, 0.1661, 0.0331 },
      { 0.0417,    0.0, 0.0438, 0.1746 },
      { 0.1744, 0.0440,    0.0, 0.0418 },
      { 0.0331, 0.1662, 0.0405,    0.0 }
    };
    model = SubstitutionModel(Q);
  }
  ~FixtureVario() {
    BOOST_TEST_MESSAGE( "teardown fixture" );
  }

  string fn_vcf;
  SubstitutionModel model;
  seqio::Genome ref_genome;
  long seed = 123456789;
  RandomNumberGenerator<> rng = RandomNumberGenerator<>(seed);
};

BOOST_FIXTURE_TEST_SUITE( vario, FixtureVario )

/* read benchmark variant set and calculate sumstats */
BOOST_AUTO_TEST_CASE( vcf_sumstats )
{
  VariantSet variant_set;
  std::map<string, vector<Genotype> > gtMatrix;
  BOOST_TEST_MESSAGE( "loading benchmark VCF file '" << fn_vcf << "'..." );
  readVcf(fn_vcf, variant_set, gtMatrix);
  vector<Variant> variants = variant_set.vec_variants;
  BOOST_TEST_MESSAGE( "  read " << variants.size() << " variants." );
  BOOST_CHECK( variants.size() == 3100839 );
  BOOST_CHECK( gtMatrix.size() == 1 );
  BOOST_CHECK( gtMatrix[0].size() == variants.size() );

  BOOST_TEST_MESSAGE( "calculating VCF sumstats..." );
  //variant_set.calculateSumstats();
  double (&f)[4][4] = variant_set.mat_freqs;

  BOOST_TEST_MESSAGE( "\tnucleotide substitution frequencies:" );
  BOOST_TEST_MESSAGE( "   |    A   |    C   |    G   |    T   " );
  BOOST_TEST_MESSAGE( format(" A | %0.4f | %0.4f | %0.4f | %0.4f ") % f[0][0] % f[0][1] % f[0][2] % f[0][3] );
  BOOST_TEST_MESSAGE( format(" C | %0.4f | %0.4f | %0.4f | %0.4f ") % f[1][0] % f[1][1] % f[1][2] % f[1][3] );
  BOOST_TEST_MESSAGE( format(" G | %0.4f | %0.4f | %0.4f | %0.4f ") % f[2][0] % f[2][1] % f[2][2] % f[2][3] );
  BOOST_TEST_MESSAGE( format(" T | %0.4f | %0.4f | %0.4f | %0.4f ") % f[3][0] % f[3][1] % f[3][2] % f[3][3] );
}

/* generate set of novel variants */
BOOST_AUTO_TEST_CASE( generate_variants )
{
  int num_vars = 1000;
  boost::timer::auto_cpu_timer t;
  BOOST_TEST_MESSAGE( "duplicating reference genome..." );
  ref_genome.duplicate();

  BOOST_TEST_MESSAGE( "generating genomic variants (using substitution frequencies)..." );
  vector<Variant> variants = generateVariants(num_vars, ref_genome, model, rng);
  BOOST_TEST_MESSAGE( "done generating 1000 variants, here are the first 10: " );
  BOOST_TEST_MESSAGE( "  id, chr, bp, ref, alt, copy" );
  for (int i=0; i<10; ++i) {
    BOOST_TEST_MESSAGE( format("%s,%s,%d,%s,%s,%s") % variants[i].id % variants[i].chr % variants[i].pos % variants[i].alleles[0] % variants[i].alleles[1] % variants[i].chr_copy );
  }

  // export variants to file
  vector<Variant> var_sorted = Variant::sortByPositionPoly(variants);
  vector<int> vec_idx = { 0 };
  vector<string> labels = { "healthy" };
  vector<vector<bool>> mm(1, vector<bool>(num_vars, true));
  auto fn_out = "ref_variants.vcf";
  ofstream fs_out;
  fs_out.open(fn_out);
  vario::writeVcf(ref_genome.records, var_sorted, vec_idx, labels, mm, fs_out);
  fs_out.close();
  BOOST_TEST_MESSAGE( " variants are in file '" << fn_out << "'." );


  BOOST_TEST_MESSAGE( "\ngenerating personal genome..." );
  BOOST_TEST_MESSAGE( "  applying variants" );
  vario::applyVariants(ref_genome, variants);
  BOOST_TEST_MESSAGE( "  writing personal genome to file 'pers.fa'" );
  ofstream fs_pg;
  fs_pg.open("pers.fa");
  writeFasta(ref_genome.records, fs_pg);
  fs_pg.close();
}

/* compare distributions for variants generated
 *  A) based on mutated nucleotide frequencies from VCF
 *  B) when positions are picked randomly
 */
BOOST_AUTO_TEST_CASE( var_dist )
{
  boost::timer::auto_cpu_timer t;
  BOOST_TEST_MESSAGE( "generating 10000 genomic variants (using substitution frequencies)..." );
  vector<Variant> vars_model = generateVariants(10000, ref_genome, model, rng);
  VariantSet varset_model;
  varset_model.vec_variants = vars_model;
  varset_model.calculateSumstats();
  double (&f)[4][4] = varset_model.mat_freqs;
  BOOST_TEST_MESSAGE( "\tnucleotide substitution frequencies:" );
  BOOST_TEST_MESSAGE( "   |    A   |    C   |    G   |    T   " );
  BOOST_TEST_MESSAGE( format(" A | %0.4f | %0.4f | %0.4f | %0.4f ") % f[0][0] % f[0][1] % f[0][2] % f[0][3] );
  BOOST_TEST_MESSAGE( format(" C | %0.4f | %0.4f | %0.4f | %0.4f ") % f[1][0] % f[1][1] % f[1][2] % f[1][3] );
  BOOST_TEST_MESSAGE( format(" G | %0.4f | %0.4f | %0.4f | %0.4f ") % f[2][0] % f[2][1] % f[2][2] % f[2][3] );
  BOOST_TEST_MESSAGE( format(" T | %0.4f | %0.4f | %0.4f | %0.4f ") % f[3][0] % f[3][1] % f[3][2] % f[3][3] );

  BOOST_TEST_MESSAGE( "generating 10000 genomic variants (random positions)..." );
  vector<Variant> vars_random = generateVariantsRandomPos(10000, ref_genome, model, rng);
  VariantSet varset_random;
  varset_random.vec_variants = vars_random;
  varset_random.calculateSumstats();
  double (&g)[4][4] = varset_random.mat_freqs;
  BOOST_TEST_MESSAGE( "\tnucleotide substitution frequencies:" );
  BOOST_TEST_MESSAGE( "   |    A   |    C   |    G   |    T   " );
  BOOST_TEST_MESSAGE( format(" A | %0.4f | %0.4f | %0.4f | %0.4f ") % g[0][0] % g[0][1] % g[0][2] % g[0][3] );
  BOOST_TEST_MESSAGE( format(" C | %0.4f | %0.4f | %0.4f | %0.4f ") % g[1][0] % g[1][1] % g[1][2] % g[1][3] );
  BOOST_TEST_MESSAGE( format(" G | %0.4f | %0.4f | %0.4f | %0.4f ") % g[2][0] % g[2][1] % g[2][2] % g[2][3] );
  BOOST_TEST_MESSAGE( format(" T | %0.4f | %0.4f | %0.4f | %0.4f ") % g[3][0] % g[3][1] % g[3][2] % g[3][3] );
}

/* test data structures for CNVs */
BOOST_AUTO_TEST_CASE( cnv )
{
  // NOTE: reference genome generated in FixtureVario()

  // build interval map for chromosomes
  // TODO: create method in seqio::Genome for this
  // imap_chrom;
}

BOOST_AUTO_TEST_SUITE_END()
