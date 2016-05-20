#include <boost/test/unit_test.hpp>

#include "random.hpp"
#include "seqio.hpp"
#include "vario.hpp"
#include <boost/format.hpp>
using boost::format;
using evolution::SubstitutionModel;
using namespace std;
using namespace vario;

struct FixtureVario {
  FixtureVario() {
    BOOST_TEST_MESSAGE( "setup fixure" );
    fn_vcf = "data/germline.vcf";
    BOOST_TEST_MESSAGE( "loading benchmark VCF file '" << fn_vcf << "'..." );
    readVcf(fn_vcf, variant_set, gtMatrix);
    BOOST_TEST_MESSAGE( "  read " << variant_set.vec_variants.size() << " variants." );

    double Q[4][4]   = {
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
  VariantSet variant_set;
  vector<vector<Genotype> > gtMatrix;
  SubstitutionModel model;
  long seed = 123456789;
};

BOOST_FIXTURE_TEST_SUITE( vario, FixtureVario )

/* read benchmark variant set */
BOOST_AUTO_TEST_CASE( vcf_sanity_check )
{
  vector<Variant> variants = variant_set.vec_variants;
  BOOST_CHECK( variants.size() == 3100839 );
  BOOST_CHECK( gtMatrix.size() == 1 );
  BOOST_CHECK( gtMatrix[0].size() == variants.size() );
}

/* calculate sumstats */
BOOST_AUTO_TEST_CASE( vcf_sumstats )
{
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
  // initialize random number generator
  RandomNumberGenerator<> rng(seed);
  BOOST_TEST_MESSAGE( "reading benchmark genome (data/Hg18.chr20-22.fa)" );
  seqio::Genome genome = seqio::Genome("data/Hg18.chr20-22.fa");
  BOOST_TEST_MESSAGE( "generating genomic variants (using substitution frequencies)..." );
  vector<Variant> variants = generateVariants(1000, genome, model, rng);
  //BOOST_TEST_MESSAGE( "done." );
}

BOOST_AUTO_TEST_SUITE_END()
