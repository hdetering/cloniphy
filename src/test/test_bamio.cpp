#include <boost/test/unit_test.hpp>
#include <boost/timer/timer.hpp>

#include "../core/clone.hpp"
#include "../core/random.hpp"
#include "../core/seqio.cpp"
using namespace seqio;
#include "../core/treeio.hpp"
using namespace treeio;
#include "../core/vario.cpp"
using namespace vario;
#include <boost/format.hpp>
using boost::format;
using boost::str;
#include <fstream>
#include <seqan/bam_io.h>
using seqan::BamAlignmentRecord;
using seqan::BamFileIn;
using seqan::BamHeader;
using seqan::CharString;
using seqan::FormattedFileContext;
using seqan::getAbsolutePath;

using namespace std;

struct FixtureBamio {
  FixtureBamio() {
    BOOST_TEST_MESSAGE( "setup fixure" );

    double Q[4][4] = {
      {    0.0, 0.0406, 0.1661, 0.0331 },
      { 0.0417,    0.0, 0.0438, 0.1746 },
      { 0.1744, 0.0440,    0.0, 0.0418 },
      { 0.0331, 0.1662, 0.0405,    0.0 }
    };
    model = SubstitutionModel(Q);
  }
  ~FixtureBamio() {
    BOOST_TEST_MESSAGE( "teardown fixture" );
  }

  long seed = 123456789;
  RandomNumberGenerator<> rng = RandomNumberGenerator<>(seed);
  SubstitutionModel model;
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
  Tree<Clone> tree(num_clones);
  tree.generateRandomTopologyLeafsOnly(random_dbl);
  tree.varyBranchLengths(random_gamma);
  int num_mutations = 1010;
  int num_transmuts =   10;
  tree.evolve(num_mutations, num_transmuts, rng);

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
  BOOST_TEST_MESSAGE( str(format("sum_probs: %0.4f") % sum_probs) );
  tree.assignWeights(w);
  // extract info about visible clones
  vector<Clone*> vec_vis_clones = tree.getVisibleNodes();
  vector<int> vec_clone_idx;
  vector<string> vec_clone_lbl;
  vector<double> vec_clone_weight;
  for (Clone *c : vec_vis_clones) {
    fprintf(stdout, "clone %d: \"%s\"\t(%.4f)\n", c->index, c->label.c_str(), c->weight);
    vec_clone_idx.push_back(c->index);
    vec_clone_lbl.push_back(c->label);
    vec_clone_weight.push_back(c->weight);
  }

  BOOST_TEST_MESSAGE( "Writing clone tree to file\n\tbulk_clone_tree.dot" );
  ofstream fs_dot;
  fs_dot.open("bulk_clone_tree.dot");
  tree.printDot(tree.m_root, fs_dot);
  fs_dot.close();

  // get mutation matrix
  BOOST_TEST_MESSAGE( "Building mutation matrix from tree..." );
  int num_nodes = tree.m_numNodes;
  vector<vector<bool>> mm(num_nodes, vector<bool>(num_mutations, false));
  tree.m_root->populateMutationMatrixRec(mm);
  // TODO: make a convenience function for this!
  BOOST_TEST_MESSAGE( "  writing matrix for visible clones to file 'mm_vis.csv'" );
  ofstream fs_mm;
  fs_mm.open("mm_vis.csv");
  for (auto cidx : vec_clone_idx) {
    fs_mm << tree.m_vecNodes[cidx]->label;
    for (auto cell : mm[cidx])
      fs_mm << "," << cell;
    fs_mm << endl;
  }
  fs_mm.close();

  // read personal genome
  Genome genome("pers.fa");
  genome.ploidy = 2;
  genome.indexRecords();

  // generate (sub)clonal variants
  BOOST_TEST_MESSAGE( "generating genomic variants (using substitution frequencies)..." );
  vector<Variant> variants = generateVariants(num_mutations, genome, model, rng);
  vector<Variant> var_sorted = Variant::sortByPositionPoly(variants);
  auto fn_out = "pers.bulk.vcf";
  ofstream fs_out;
  fs_out.open(fn_out);
  vector<int> vec_idx(vec_vis_clones.size());
  iota(vec_idx.begin(), vec_idx.end(), 0);
  vario::writeVcf(genome.records, var_sorted, vec_idx, vec_clone_lbl, mm, fs_out);
  fs_out.close();
  BOOST_TEST_MESSAGE( " variants are in file '" << fn_out << "'." );

  /*------------------*
   * process BAM file *
   *------------------*/

  function<short()> rand_copy = rng.getRandomFunctionInt(short(0), short(genome.ploidy-1));
  function<int()> rand_idx = rng.getRandomIndexWeighted(vec_clone_weight);
  auto it_var_mark = var_sorted.begin(); // the leftmost variant we need to consider
  CharString bamFileName = getAbsolutePath("build/pers_reads.sam");
  BamFileIn bamFileIn;
  if (!open(bamFileIn, toCString(bamFileName))) {
    std::cerr << "ERROR: Could not open " << bamFileName << std::endl;
  };
  // NOTE: reading the header first is MANDATORY!
  BamHeader header;
  readHeader(header, bamFileIn);
  // just for fun: look at the header content
  typedef FormattedFileContext<BamFileIn, void>::Type TBamContext;
  TBamContext const & bamContext = context(bamFileIn);
  for (unsigned i = 0; i < length(contigNames(bamContext)); ++i)
    std::cout << contigNames(bamContext)[i] << '\t'
              << contigLengths(bamContext)[i] << '\n';

  ofstream fs_idx, fs_fq, fs_log;
  fs_idx.open("idx.csv");
  fs_fq.open("bulk_reads.fq");
  fs_log.open("bamio_bulk.log");
  BamAlignmentRecord read1, read2;
  while (!atEnd(bamFileIn)) {
    readRecord(read1, bamFileIn);
    readRecord(read2, bamFileIn);
    // pick random clone (weighted) to which next read pair belongs
    int c_idx = vec_clone_idx[rand_idx()];
    int r1_begin = read1.beginPos;
    int r2_begin = read2.beginPos;
    int r1_len = getAlignmentLengthInRef(read1);
    int r2_len = getAlignmentLengthInRef(read2);
    CharString r1_ref = contigNames(bamContext)[read1.rID];
    CharString r2_ref = contigNames(bamContext)[read2.rID];
    char r1_rc = seqan::hasFlagRC(read1) ? '+' : '-';
    char r2_rc = seqan::hasFlagRC(read2) ? '+' : '-';
    string r1_qual = toCString(read1.qual);
    string r2_qual = toCString(read2.qual);
    fs_idx << str(format("%s:%s(%d,%d) -- %s:%s(%d,%d)\n")
      % r1_ref % r1_rc % r1_begin % (r1_begin+r1_len)
      % r2_ref % r2_rc % r2_begin % (r2_begin+r2_len));

    // modify read pair to match assigned clone
    read1.qName += str(format("-%d/1") % c_idx);
    read2.qName += str(format("-%d/2") % c_idx);

    auto it_var = var_sorted.begin();
    // advance variant iterator to position in current chromosome
    while (it_var->chr != toCString(r1_ref))
      ++it_var;
    // advance variant iterator to first position past begin of first read
    while (it_var->pos < r1_begin && (it_var+1)->chr == toCString(r1_ref))
      ++it_var;
    // identify variables affecting read pair
    while (it_var->pos < r2_begin+r2_len && it_var->chr == toCString(r1_ref) && it_var != var_sorted.end()) {
      if (!mm[c_idx][it_var->idx_mutation]) { // does read pair's clone carry mutation?
        ++it_var;
        continue;
      }
      int r1_var_pos = it_var->pos - r1_begin;
      if (r1_var_pos >= 0 && r1_var_pos < r1_len) { // mutate read1
        fs_log << str(format("%s:%d\t%s->%s\n") % toCString(read1.qName) % r1_var_pos % it_var->alleles[0].c_str() % it_var->alleles[1].c_str());
        read1.seq[r1_var_pos] = it_var->alleles[1][0];
      }
      int r2_var_pos = it_var->pos - r2_begin;
      if (r2_var_pos >= 0 && r2_var_pos < r2_len) { // mutate read2
        fs_log << str(format("%s:%d\t%s->%s\n") % toCString(read2.qName) % r2_var_pos % it_var->alleles[0].c_str() % it_var->alleles[1].c_str());
        read2.seq[r2_var_pos] = it_var->alleles[1][0];
      }
      ++it_var;
    }

    // re-reverse complement reads on opposite strand (to avoid strand bias)
    if (seqan::hasFlagRC(read1)) {
      seqan::reverseComplement(read1.seq);
      reverse(r1_qual.begin(), r1_qual.end());
    }
    if (seqan::hasFlagRC(read2)) {
      seqan::reverseComplement(read2.seq);
      reverse(r2_qual.begin(), r2_qual.end());
    }

    // write (potentially modified) reads to FASTQ file
    fs_fq << str(format("@%s\n%s\n+\n%s\n") % read1.qName % read1.seq % r1_qual);
    fs_fq << str(format("@%s\n%s\n+\n%s\n") % read2.qName % read2.seq % r2_qual);

    //fs_idx << r1_begin << " (" << r1_len << ")" <<  "  " << r2_begin << " (" << r2_len << ")" << endl;
    //fprintf(stderr, "%d\n", c_idx);
  }
  fs_idx.close();
  fs_fq.close();
  fs_log.close();


  /*auto fn_reads1 = "reads1.fq";
  BOOST_TEST_MESSAGE( "Writing results to output file" );
  BOOST_TEST_MESSAGE( str(format("  %s") % fn) );
  ofstream fs;
  fs.open(fn);
  for (double v : vec_values) {
    fs << format("%f") % v << endl;
  }
  fs.close();*/
}

BOOST_AUTO_TEST_SUITE_END()
