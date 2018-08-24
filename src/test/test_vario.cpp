#include <boost/test/unit_test.hpp>
#include <boost/timer/timer.hpp>

#include "../core/random.hpp"
#include "../core/seqio.hpp"
#include "../core/vario.hpp"
#include "../core/vario/VariantStore.hpp"
#include <boost/icl/interval_map.hpp>
using namespace boost::icl;
#include <fstream>
#include <vector>
using stringio::format;
using evolution::GermlineSubstitutionModel;
using namespace std;
using namespace seqio;
using namespace vario;

struct FixtureVario {
  FixtureVario() {
    BOOST_TEST_MESSAGE( "setup fixure" );
    fn_vcf = "data/germline.vcf";

    unsigned long ref_genome_len = 1000000;
    // int num_ref_seq = 10;
    // int seq_len_mean = 100000;
    // int seq_len_sd = 10000;
    int num_ref_seq = 1;
    int seq_len_mean = 100000;
    int seq_len_sd = 0;
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
    model = GermlineSubstitutionModel(Q);
  }
  ~FixtureVario() {
    BOOST_TEST_MESSAGE( "teardown fixture" );
  }

  string fn_vcf;
  GermlineSubstitutionModel model;
  seqio::GenomeReference ref_genome;
  long seed = 123456789;
  RandomNumberGenerator rng = RandomNumberGenerator(seed);
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
  BOOST_TEST_MESSAGE( format(" A | %0.4f | %0.4f | %0.4f | %0.4f ", f[0][0], f[0][1], f[0][2], f[0][3]) );
  BOOST_TEST_MESSAGE( format(" C | %0.4f | %0.4f | %0.4f | %0.4f ", f[1][0], f[1][1], f[1][2], f[1][3]) );
  BOOST_TEST_MESSAGE( format(" G | %0.4f | %0.4f | %0.4f | %0.4f ", f[2][0], f[2][1], f[2][2], f[2][3]) );
  BOOST_TEST_MESSAGE( format(" T | %0.4f | %0.4f | %0.4f | %0.4f ", f[3][0], f[3][1], f[3][2], f[3][3]) );
}

/* generate set of novel variants */
/* TODO: Test has to be rewritten (use VariantStore). */
BOOST_AUTO_TEST_CASE( germline )
{
  int num_vars = 1000;
  boost::timer::auto_cpu_timer t;
  VariantStore var_store;

  BOOST_TEST_MESSAGE( "generating genomic variants (using substitution frequencies)..." );
  var_store.generateGermlineVariants(num_vars, ref_genome, model, 0.1, rng);
  vector<Variant> variants = var_store.getGermlineSnvVector();

  BOOST_TEST_MESSAGE( "done generating 1000 variants, here are the first 10: " );
  BOOST_TEST_MESSAGE( "  id, chr, bp, ref, alt" );
  for (int i=0; i<10; ++i) {
    BOOST_TEST_MESSAGE( format("%s,%s,%d,%s,%s", variants[i].id, variants[i].chr, variants[i].pos, variants[i].alleles[0], variants[i].alleles[1]) );
  }

  // export variants to file
  vector<Variant> var_sorted = Variant::sortByPositionLex(variants);
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
  // TODO: use VariantStore object!
  //vario::applyVariants(ref_genome, variants);
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
  VariantStore var_store;

  double mut_gl_hom = 0.1; // ratio of homozygous variants

  BOOST_TEST_MESSAGE( "generating 10000 genomic variants (using substitution frequencies)..." );
  var_store.generateGermlineVariants(10000, ref_genome, model, mut_gl_hom, rng);
  vector<Variant> vars_model = var_store.getGermlineSnvVector();
  VariantSet varset_model;
  varset_model.vec_variants = vars_model;
  varset_model.calculateSumstats();
  double (&f)[4][4] = varset_model.mat_freqs;
  BOOST_TEST_MESSAGE( "\tnucleotide substitution frequencies:" );
  BOOST_TEST_MESSAGE( "   |    A   |    C   |    G   |    T   " );
  BOOST_TEST_MESSAGE( format(" A | %0.4f | %0.4f | %0.4f | %0.4f ", f[0][0], f[0][1], f[0][2], f[0][3]) );
  BOOST_TEST_MESSAGE( format(" C | %0.4f | %0.4f | %0.4f | %0.4f ", f[1][0], f[1][1], f[1][2], f[1][3]) );
  BOOST_TEST_MESSAGE( format(" G | %0.4f | %0.4f | %0.4f | %0.4f ", f[2][0], f[2][1], f[2][2], f[2][3]) );
  BOOST_TEST_MESSAGE( format(" T | %0.4f | %0.4f | %0.4f | %0.4f ", f[3][0], f[3][1], f[3][2], f[3][3]) );

  BOOST_TEST_MESSAGE( "generating 10000 genomic variants (random positions)..." );
  vector<Variant> vars_random = generateVariantsRandomPos(10000, ref_genome, model, rng);
  VariantSet varset_random;
  varset_random.vec_variants = vars_random;
  varset_random.calculateSumstats();
  double (&g)[4][4] = varset_random.mat_freqs;
  BOOST_TEST_MESSAGE( "\tnucleotide substitution frequencies:" );
  BOOST_TEST_MESSAGE( "   |    A   |    C   |    G   |    T   " );
  BOOST_TEST_MESSAGE( format(" A | %0.4f | %0.4f | %0.4f | %0.4f ", g[0][0], g[0][1], g[0][2], g[0][3]) );
  BOOST_TEST_MESSAGE( format(" C | %0.4f | %0.4f | %0.4f | %0.4f ", g[1][0], g[1][1], g[1][2], g[1][3]) );
  BOOST_TEST_MESSAGE( format(" G | %0.4f | %0.4f | %0.4f | %0.4f ", g[2][0], g[2][1], g[2][2], g[2][3]) );
  BOOST_TEST_MESSAGE( format(" T | %0.4f | %0.4f | %0.4f | %0.4f ", g[3][0], g[3][1], g[3][2], g[3][3]) );
}

/* test data structures for CNVs */
BOOST_AUTO_TEST_CASE( cnv )
{
  // NOTE: reference genome generated in FixtureVario()
  unsigned n_ref_chr = ref_genome.num_records;
  GenomeInstance g_inst(ref_genome);
  BOOST_TEST_MESSAGE( "Number of chromosomes in GenomeInstance: " << g_inst.map_id_chr.size() );
  // check: does genome instance have same number of chromomsomes as reference?
  BOOST_CHECK( g_inst.map_id_chr.size() == n_ref_chr );
  BOOST_CHECK( g_inst.vec_chr.size() == 2*n_ref_chr );
  BOOST_CHECK( g_inst.vec_chr_len.size() == 2*n_ref_chr );

  VariantStore var_store;

  // Initialize mutations
  Mutation mut_cnv_del;
  mut_cnv_del.id = 0;
  mut_cnv_del.is_cnv = true;
  CopyNumberVariant var_cnv_del;
  var_cnv_del.is_deletion = true;
  var_cnv_del.is_telomeric = false;
  var_cnv_del.is_forward = true;
  var_cnv_del.ref_chr = "chr0";
  var_cnv_del.len_rel = 0.5;
  var_cnv_del.start_rel = 0.25;

  Mutation mut_cnv_amp;
  mut_cnv_amp.id = 1;
  mut_cnv_amp.is_cnv = true;
  CopyNumberVariant var_cnv_amp;
  var_cnv_amp.is_deletion = false;
  var_cnv_amp.is_telomeric = false;
  var_cnv_amp.is_forward = true;
  var_cnv_amp.ref_chr = "chr0";
  var_cnv_amp.len_rel = 0.2;
  var_cnv_amp.start_rel = 0.6;

  Mutation mut_cnv_wgd;
  mut_cnv_wgd.id = 2;
  mut_cnv_wgd.is_cnv = true;
  CopyNumberVariant var_cnv_wgd;
  var_cnv_wgd.is_wgd = true;

  Mutation mut_snv_1;
  mut_snv_1.id = 3;
  mut_snv_1.is_snv = true;
  Variant var_snv_1("snv1", "chr0", 100);

  // register variants in VariantStore
  var_store.map_id_cnv[mut_cnv_del.id] = var_cnv_del;
  var_store.map_id_cnv[mut_cnv_amp.id] = var_cnv_amp;
  var_store.map_id_cnv[mut_cnv_wgd.id] = var_cnv_wgd;
  var_store.map_id_snv[mut_snv_1.id] = var_snv_1;
  shared_ptr<ChromosomeInstance> chr_src = g_inst.map_id_chr[var_snv_1.chr][0];
  SegmentCopy seg_src = chr_src->lst_segments.front();
  var_store.map_seg_vars[seg_src.id] = { mut_snv_1.id };

  BOOST_TEST_MESSAGE( "Genome initially:\n" << g_inst );

  // FOCAL DELETION
  //-------------------------------------

  BOOST_TEST_MESSAGE( "\nNow performing FOCAL DELETION..." );
  BOOST_TEST_MESSAGE( "--------------------------------\n" );

  BOOST_TEST_MESSAGE( "CNV: DEL<chr=" << var_cnv_del.ref_chr << ", start=" << var_cnv_del.start_rel << ", len=" << var_cnv_del.len_rel << ">" );

  var_store.applyMutation(mut_cnv_del, g_inst, rng);
  BOOST_TEST_MESSAGE( "Genome after focal deletion:\n" << g_inst );

  // FOCAL AMPLIFICATION
  //-------------------------------------

  BOOST_TEST_MESSAGE( "\nNow performing FOCAL AMPLIFICATION..." );
  BOOST_TEST_MESSAGE( "-------------------------------------\n" );

  BOOST_TEST_MESSAGE( "CNV: CPY<chr=" << var_cnv_amp.ref_chr << ", start=" << var_cnv_amp.start_rel << ", len=" << var_cnv_amp.len_rel << ">" );

  var_store.applyMutation(mut_cnv_amp, g_inst, rng);
  BOOST_TEST_MESSAGE( "Genome after focal amplification:\n" << g_inst );


  // WHOLE GENOME DUPLICATION (WGD)
  //-------------------------------------

  BOOST_TEST_MESSAGE( "\nNow performing WHOLE GENOME DUPLICATION..." );

  var_store.applyMutation(mut_cnv_wgd, g_inst, rng);
  BOOST_TEST_MESSAGE( "Genome after WGD:\n" << g_inst );

  // where chromosomes duplicated?
  BOOST_CHECK( g_inst.map_id_chr.size() == n_ref_chr );
  BOOST_CHECK( g_inst.vec_chr.size() == 4*n_ref_chr );
  BOOST_CHECK( g_inst.vec_chr_len.size() == 4*n_ref_chr );

  // were variants duplicated correctly?
  shared_ptr<ChromosomeInstance> chr_copy = g_inst.map_id_chr[var_snv_1.chr][2];
  SegmentCopy seg_copy = chr_copy->lst_segments.front();
  BOOST_CHECK( var_store.map_seg_vars.count(seg_copy.id) == 1 );
  auto v_mut_copy = var_store.map_seg_vars[seg_copy.id];
  BOOST_CHECK( find(v_mut_copy.begin(), v_mut_copy.end(), mut_snv_1.id) != v_mut_copy.end() );
}

BOOST_AUTO_TEST_SUITE_END()
