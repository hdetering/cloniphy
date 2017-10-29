#include <boost/test/unit_test.hpp>
#include <boost/timer/timer.hpp>

#include "../core/bamio.hpp"
#include "../bamio/BulkSampleGenerator.hpp"
using namespace vario;
#include "../core/clone.hpp"
#include "../core/config/ConfigStore.hpp"
#include "../core/random.hpp"
#include "../core/seqio.cpp"
using namespace seqio;
#include "../core/treeio.hpp"
using namespace treeio;
#include <boost/format.hpp>
using boost::str;
#include <fstream>
#include <sstream>

using namespace std;

struct FixtureBamio {
  FixtureBamio() {
    BOOST_TEST_MESSAGE( "setup fixure" );

    // init germline mutation model
    double Q[4][4] = {
      {    0.0, 0.0406, 0.1661, 0.0331 },
      { 0.0417,    0.0, 0.0438, 0.1746 },
      { 0.1744, 0.0440,    0.0, 0.0418 },
      { 0.0331, 0.1662, 0.0405,    0.0 }
    };
    model_gl = GermlineSubstitutionModel(Q);

    // init somatic mutation model
    string fn_somatic = "../resources/signatures_probabilities.txt";
    map<string, double> contrib;
    contrib["Signature 1"] = 0.2;
    contrib["Signature 22"] = 0.8;
    model_sm = SomaticSubstitutionModel(fn_somatic, contrib);

    char* argv[3] = {const_cast<char*>("test_bamio"), const_cast<char*>("-c"), const_cast<char*>("../config.yml")};
    config.parseArgs(3, argv);
  }
  ~FixtureBamio() {
    BOOST_TEST_MESSAGE( "teardown fixture" );
  }

  long seed = 123456789;
  RandomNumberGenerator<> rng = RandomNumberGenerator<>(seed);
  GermlineSubstitutionModel model_gl;
  SomaticSubstitutionModel model_sm;
  SomaticCnvModel model_cnv;
  VariantStore var_store;
  config::ConfigStore config;
};

BOOST_FIXTURE_TEST_SUITE( bamio, FixtureBamio )

/* generate random gamma-distributed values */
BOOST_AUTO_TEST_CASE( bulk )
{
  // load genealogy of 5 clones
  /*auto fn_tree = "random_sample_tree.tre";
  node root;
  readNewick(fn_tree, root);
  Tree<Clone> tree(root);*/

  // TODO: this should be replaced by loading a tree from file (see above)
  // -> integrate flags "is_visible" and "weight" into Newick format
  function<double()> random_dbl = rng.getRandomFunctionDouble(0.0, 1.0);
  function<double()> random_gamma = rng.getRandomGamma(2.0, 0.25);
  int num_clones = 5;
  treeio::Tree<Clone> tree(num_clones);
  tree.generateRandomTopologyLeafsOnly(random_dbl);
  tree.varyBranchLengths(random_gamma);
  int num_mutations = 1010;
  int num_transmuts =   10;
  tree.dropSomaticMutations(num_mutations, num_transmuts, rng);
  // TODO: assignMutationType

  BOOST_CHECK( tree.m_numVisibleNodes == num_clones+1 );
  BOOST_CHECK( tree.m_numNodes == tree.m_vecNodes.size() );

  // assign random weights
  /*vector<double> w = rng.getRandomProbs(num_clones);
  double c = 0.1; // contamination with normal cells
  for (auto &p : w) p*=(1-c);
  w.insert(w.begin(), c);*/

  // assign user-defined weights (cellular prevalence)
  vector<double> w = { 0.1, 0.4, 0.3, 0.1, 0.075, 0.025 };
  double sum_probs = 0.0;
  for (auto p : w) sum_probs += p;
  BOOST_CHECK( sum_probs == 1.0 );
  BOOST_TEST_MESSAGE( str(boost::format("sum_probs: %0.4f") % sum_probs) );
  tree.assignWeights(w);

  BOOST_TEST_MESSAGE( "Writing clone tree to file\n\tbulk_clone_tree.dot" );
  tree.printDot("bulk_clone_tree.dot");
  BOOST_TEST_MESSAGE( "Writing clone tree to file\n\tbulk_clone_tree.nwk" );
  tree.printNewick("bulk_clone_tree.nwk");
  return;

  // get mutation matrix
  BOOST_TEST_MESSAGE( "Building mutation matrix from tree..." );
  int num_nodes = tree.m_numNodes;
  vector<vector<bool>> mm(num_nodes, vector<bool>(num_mutations, false));
  tree.m_root->populateMutationMatrixRec(mm);
  // write mutation matrix to file
  string fn_mm = "mm_vis.csv";
  BOOST_TEST_MESSAGE( "  writing matrix for visible clones to file '" << fn_mm << "'" );
  tree.writeMutationMatrix(fn_mm);

  // read personal genome
  GenomeReference genome("pers.fa");
  genome.ploidy = 2;
  genome.indexRecords();

  // generate (sub)clonal variants
  BOOST_TEST_MESSAGE( "generating genomic variants (using substitution frequencies)..." );
  vector<Mutation> vec_mutations = vector<Mutation>(num_mutations);
  for (auto &m : vec_mutations) {
    m.is_snv = true;
    m.is_cnv = false;
  }
  var_store.generateSomaticVariants(vec_mutations, genome, model_sm, model_cnv, rng);
  vector<Variant> variants = var_store.getSomaticSnvVector();
  vector<Variant> var_sorted = Variant::sortByPositionLex(variants);
  auto fn_out = "pers.bulk.vcf";
  ofstream fs_out;
  fs_out.open(fn_out);
  vector<shared_ptr<Clone>> vec_vis_clones = tree.getVisibleNodes();
  vector<int> vec_idx;
  vector<string> vec_lbl;
  for (auto clone : vec_vis_clones) {
    vec_idx.push_back(clone->index);
    vec_lbl.push_back(clone->label);
  }

  vario::writeVcf(genome.records, var_sorted, vec_idx, vec_lbl, mm, fs_out);
  fs_out.close();
  BOOST_TEST_MESSAGE( " variants are in file '" << fn_out << "'." );

  /* process BAM file */
  VariantSet varset(variants);
  bamio::mutateReads("bulk_reads.fq", "bulk_reads.sam", "build/pers_reads.sam", varset, tree, w, "BULKSIM", 2, rng);
}

BOOST_AUTO_TEST_CASE( multisample )
{
  // get parameters from config file
  int num_clones = config.getValue<int>("clones");
  map<string, vector<double>> sample_mtx = config.getMatrix<double>("samples");
  int num_mutations = config.getValue<int>("mutations");
  int num_transmuts = config.getValue<int>("init-muts");

  // initialize random functions
  function<double()> random_dbl = rng.getRandomFunctionDouble(0.0, 1.0);
  function<double()> random_gamma = rng.getRandomGamma(2.0, 0.25);

  // generate clone tree
  treeio::Tree<Clone> tree(num_clones);
  tree.generateRandomTopologyLeafsOnly(random_dbl);
  tree.varyBranchLengths(random_gamma);
  // drop mutations on tree
  tree.dropSomaticMutations(num_mutations, num_transmuts, rng);
  // TODO: assignMutationType
  BOOST_CHECK( tree.m_numVisibleNodes == num_clones+1 );
  BOOST_CHECK( tree.m_numNodes == tree.m_vecNodes.size() );
  BOOST_TEST_MESSAGE( "Writing clone tree to file\n\tmultisample.tree.dot" );
  tree.printDot("multisample.tree.dot");
  BOOST_TEST_MESSAGE( "Writing clone tree to file\n\tmultisample.tree.nwk" );
  tree.printNewick("multisample.tree.nwk");

  // get mutation matrix
  BOOST_TEST_MESSAGE( "Building mutation matrix from tree..." );
  int num_nodes = tree.m_numNodes;
  vector<vector<bool>> mm(num_nodes, vector<bool>(num_mutations, false));
  tree.m_root->populateMutationMatrixRec(mm);
  // write mutation matrix
  string fn_mm = "multisample.mm.csv";
  BOOST_TEST_MESSAGE( "Writing mutation matrix for visible clones to file '" << fn_mm << "'" );
  tree.writeMutationMatrix(fn_mm);

  // read personal genome
  GenomeReference genome("pers.fa");
  genome.ploidy = 2;
  genome.indexRecords();

  // generate (sub)clonal variants
  BOOST_TEST_MESSAGE( "generating somatic variants (using substitution frequencies)..." );
  vector<Mutation> vec_mutations = vector<Mutation>(num_mutations);
  for (auto &m : vec_mutations) {
    m.is_snv = true;
    m.is_cnv = false;
  }
  var_store.generateSomaticVariants(vec_mutations, genome, model_sm, model_cnv, rng);
  vector<Variant> variants = var_store.getSomaticSnvVector();
  vector<Variant> var_sorted = Variant::sortByPositionLex(variants);
  vector<shared_ptr<Clone>> vec_vis_clones = tree.getVisibleNodes();
  vector<int> vec_idx;
  vector<string> vec_lbl;
  for (auto clone : vec_vis_clones) {
    vec_idx.push_back(clone->index);
    vec_lbl.push_back(clone->label);
  }

  auto fn_out = "multisample.bulk.vcf";
  ofstream fs_out;
  fs_out.open(fn_out);
  vario::writeVcf(genome.records, var_sorted, vec_idx, vec_lbl, mm, fs_out);
  fs_out.close();
  BOOST_TEST_MESSAGE( " variants are in file '" << fn_out << "'." );

  // create sequencing reads for samples by applying variants to "germline" reads
  VariantSet varset(variants);
  for (auto row_sample : sample_mtx) {
    string lbl_sample = row_sample.first;
    vector<double> w = row_sample.second;
    string fn_fastq = str(boost::format("%s.fq") % lbl_sample);
    string fn_sam = str(boost::format("%s.sam") % lbl_sample);
    // TODO: use BulkSampleGenerator object!
    //bamio::mutateReads(fn_fastq, fn_sam, "build/pers.sam", varset, tree, w, lbl_sample, genome.ploidy, rng);
  }
  BOOST_TEST_MESSAGE( "EOT" );
}

BOOST_AUTO_TEST_CASE( crc )
{
  // read clone tree from file
  treeio::Tree<Clone> tree("data/crc/clones.crc.tree");
  vector<shared_ptr<Clone>> vec_clones = tree.getVisibleNodes();

  // read variants from VCF
  string fn_vcf = "crc.somatic.vcf";
  vario::VariantSet varset;
  map<string, vector<vario::Genotype> > mtx_gt;
  vario::readVcf(fn_vcf, varset, mtx_gt);

  // read mutation map from file
  unsigned n_rows_mm = 0;
  map<string, vector<bool>> mm;
  n_rows_mm = vario::readMutMapFromCSV(mm, "crc.mm.csv");

  // process BAM file
  string fn_fq_out = "";
  string fn_sam_out = "buid/crc.test.sam";
  string fn_sam_in = "build/crc.ref.sam";
  vector<double> weights = { 0.50, 0.45, 0.05, 0.0, 0.0, 0.0, 0.0 };
  string id_sample = "T1";
  short ploidy = 2;
  bool do_write_fastq = false;

  mutateReads(fn_fq_out, fn_sam_out, fn_sam_in, varset, vec_clones, mm,
              weights, id_sample, ploidy, rng, do_write_fastq);
}

/* Test BulkSampleGenerator */
BOOST_AUTO_TEST_CASE ( bsg )
{
  using boost::filesystem::path;
  path workdir = "/home/harry/code/cloniphy/build/cnvsim";
  path fn_ref = workdir/"ref.fa";
  path path_bam = workdir/"bam";
  
  // 1. read reference genome
  BOOST_TEST_MESSAGE( "\nReading reference genome from file:\n" + fn_ref.string() );
  GenomeReference ref_genome(fn_ref.c_str());
  BOOST_TEST_MESSAGE( "Read " << ref_genome.num_records << "records:\n" );
  BOOST_TEST_MESSAGE( ref_genome );

  // 2. generate germline variants
  vario::VariantStore var_store;
  

  // 2. initialize BulkSampleGenerator
  BOOST_TEST_MESSAGE( "\nInitializing BulkSampleGenerator.\n" );
  BulkSampleGenerator bulk_gen;
  bulk_gen.initializeRefSeqs(ref_genome);

  // 3. generate read groups to assign read pairs to clones
  BOOST_TEST_MESSAGE( "\nGenerating Read Groups.\n" );
  vector<seqan::BamHeaderRecord> vec_rg;
  vector<string> vec_tag_id = {"A","B","C","D"};
  bulk_gen.generateReadGroups(vec_rg, "sample1", vec_tag_id, "Illumina", "HiSeq2500");

  // 4. merge existing reads (won't work without BAM files in path_bam)
  BOOST_TEST_MESSAGE( "\nNow merging BAM files...\n" );
  bulk_gen.mergeBulkSeqReads(path_bam, "sample1", vec_rg, var_store, rng);
}

/* Test BOOST interval container library */
BOOST_AUTO_TEST_CASE ( icl )
{
  typedef set<string> SegmentsT;
  typedef boost::icl::interval_map<int, SegmentsT> SegMapT;
  typedef boost::icl::interval_set<int> SegSetT;

  // interval map will store sequence segments [start, end)
  SegMapT imap;
  // add intervals [0,10), [0,4), [2,6)
  imap += make_pair(interval<int>::right_open(0,10), SegmentsT({"seg.1"}));
  imap += make_pair(interval<int>::right_open(0,4), SegmentsT({"seg.1.1"}));
  imap += make_pair(interval<int>::right_open(2,6), SegmentsT({"seg.1.2"}));

  cout << "Interval map content:" << endl;
  cout << imap << endl << endl;

  // look for intervals in the map which overlap different interval sets

  // first set: should overlap with {"seg.1"}
  SegSetT iset1;
  iset1.add(interval<int>::right_open(7,8));
  SegMapT isect1 = imap & iset1;
  cout << "Overlap with [7,8):" << endl;
  cout << isect1 << endl;

  // second set: should overlap with {"seg.1", "seg.1.1", "seg.1.2"}
  SegSetT iset2;
  iset2.add(interval<int>::right_open(3,4));
  SegMapT isect2 = imap & iset2;
  cout << "Overlap with [3,4):" << endl;
  cout << isect2 << endl;

  // third set: should overlap with {"seg.1", "seg.1.1"}
  SegSetT iset3;
  iset3.add(interval<int>::right_open(1,2));
  SegMapT isect3 = imap & iset3;
  cout << "Overlap with [1,2):" << endl;
  cout << isect3 << endl;
}

BOOST_AUTO_TEST_SUITE_END()
