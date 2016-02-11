/**
 * Simulation of clonal DNA sequences based on the process of somatic evolution.
 *
 * @author: Harald Detering (harald.detering@gmail.com)
 *
 */
#include "clone.hpp"
#include "seqio.hpp"
#include "treeio.hpp"
#include "vario.hpp"
#include "basicclonetree.hpp"
#include "coalescentclonetree.hpp"

#include <boost/bind.hpp>
#include <boost/format.hpp>
#include <boost/function.hpp>
#include <boost/program_options.hpp>
#include <boost/random.hpp>
// Choosing the random number generator. (mt19937: Mersenne-Twister)
typedef boost::mt19937 base_generator_type;
#include <ctime>
#include <exception>
#include <fstream>
#include <map>
#include <math.h>
#include <sstream>
#include <string>
#include <sys/stat.h>

#define PROGRAM_NAME "CloniPhy 0.1"

using namespace std;
using seqio::SeqRecord;
using seqio::Genome;
using vario::Genotype;
using vario::Variant;
using evolution::SubstitutionModel;

boost::function<float()> initRandomNumberGenerator(long seed);
bool parseArgs (int ac, char* av[], int& n_clones, std::vector<float>& freqs, int& n_mut, int& n_transmut, string& ref, string& ref_vcf, string& tree, bool verbose=true);
bool fileExists(string filename);

int main (int argc, char* argv[])
{
  // params specified by command line
  int num_clones = 0;
  int num_mutations = 0;
  int num_transmuts = 0;
  vector<float> freqs;
  string reference = "";
  string ref_vcf = "";
  string tree_fn = "";
  long seed;
  double titv = 0.5; // TODO: make this a user parameter (?)
  float ado_pct = 0.1; // TODO: make this a user parameter
  int ado_frag_len = 1000; // TODO: make this a user parameter

  // internal vars
  map<Clone*, string> clone2fn; // stores genome file names for each clone

  bool args_ok = parseArgs(argc, argv, num_clones, freqs, num_mutations, num_transmuts, reference, ref_vcf, tree_fn);
  if (!args_ok) { return EXIT_FAILURE; }

  // specify random seed
  seed = 123456789;
  //seed = time(NULL) + clock();
#ifdef DEBUG
    fprintf(stderr, "seed: %ld\n", seed);
#endif
  boost::function<float()> random = initRandomNumberGenerator(seed);
  // take that baby for a spin
  //for (int i=0; i<10; i++) { fprintf(stderr, "%.10f\n", random()); }

  Tree<Clone> tree;
  if (tree_fn.size()>0) {
    std::cerr << "Reading tree from file '" << tree_fn << "'" << std::endl;
    treeio::node root_node;
    treeio::readNewick(tree_fn, root_node);
    tree = *(new Tree<Clone>(root_node));
    num_clones = tree.getVisibleNodes().size();
    std::cerr << "num_nodes:\t" << tree.m_numNodes << std::endl;
    std::cerr << "num_clones:\t" << num_clones << std::endl;

    // if number of mutations have not been supplied specifically,
    // branch lengths are interpreted as number of mutations
    if (num_mutations == 0) {
      fprintf(stderr, "\nNumber of mutations not specified, using tree branch lengths (expected number of mutations).\n");
      // TODO: should absolute number of mutations be encodable in branch lengths?
      double dbl_branch_length = tree.getTotalBranchLength();
      num_mutations = floor(dbl_branch_length);
      fprintf(stderr, "Simulating a total of %d mutations.\n", num_mutations);
    }
  }
  else {
    tree = Tree<Clone>(num_clones);
    fprintf(stderr, "\nGenerating random topology...\n");
    tree.generateRandomTopology(random);
    fprintf(stderr, "done.\n");
  }
  fprintf(stderr, "\nNewick representation of underlying tree:\n");
  CloneTree::printNewick(tree.m_root, cerr);
  fprintf(stderr, "\n");

  tree.evolve(num_mutations, num_transmuts, random);
  fprintf(stderr, "\nNewick representation of mutated tree:\n");
  CloneTree::printNewick(tree.m_root, cerr);
  fprintf(stderr, "\n");

  fprintf(stderr, "Writing mutated tree to file 'clone_tree_mutated.dot'.\n");
  ofstream dotFileMut;
  dotFileMut.open("clone_tree_mutated.dot");
  tree.printDot(tree.m_root, dotFileMut);
  dotFileMut.close();

  /*tree.collapseZeroBranches(tree.getRoot());
  fprintf(stderr, "\nNewick representation of collapsed tree:\n");
  CoalescentCloneTree::printNewick(tree.getRoot(), cerr);
  fprintf(stderr, "\n");

  fprintf(stderr, "Writing collapsed tree to file 'clone_tree_collapsed.dot'.\n");
  ofstream dotFileCol;
  dotFileCol.open("clone_tree_collapsed.dot");
  tree.printDot(tree.getRoot(), dotFileCol);
  dotFileCol.close();*/

  // read reference sequence
  fprintf(stderr, "\nReading reference from file '%s'...", reference.c_str());
  Genome ref_genome = Genome(reference.c_str());
  fprintf(stderr, "read (%u bp in %u sequences).\n", ref_genome.length, ref_genome.num_records);
  // duplicate genome (all loci homozygous reference)
  ref_genome.duplicate();

  // if a VCF file was provided, apply germline variants
  if (ref_vcf.size() > 0) {
    fprintf(stderr, "applying germline   variants (from %s).\n", ref_vcf.c_str());
    vector<Variant> ref_variants;
    vector<vector<Genotype > > ref_gt_matrix;
    vario::readVcf(ref_vcf, ref_variants, ref_gt_matrix);
    vario::applyVariants(ref_genome, ref_variants, ref_gt_matrix[0]);
  }

  // write "healthy" genome to file
  ofstream f_fasta;
  f_fasta.open("healthy_genome.fa");
  seqio::writeFasta(ref_genome.records, f_fasta);
  f_fasta.close();
  clone2fn[tree.m_root] = "healthy_genome.fa";

  //fprintf(stderr, "generating FASTA index.\n");
  //SeqIO::indexFasta(reference.c_str());

  vector<Mutation> mutations;
  if (num_mutations > 0) {
    // generate point mutations (relative position + chr copy)
    mutations = vario::generateMutations(num_mutations, random);
    fprintf(stderr, "\nTotal set of mutations (id, rel_pos, copy):\n");
    for (int i=0; i<num_mutations; i++)
      fprintf(stderr, "%d\t%f\t%d\n", mutations[i].id, mutations[i].relPos, mutations[i].copy);
  }
  // apply mutations, creating variants in the process
  vector<Variant> variants = vector<Variant>();
  SubstitutionModel model = SubstitutionModel(ref_genome.nuc_freq, titv);
  vario::applyMutations(mutations, ref_genome, model, random, variants);

  // initialize mutation matrix
  vector<vector<short> > mutMatrix(tree.m_numNodes, std::vector<short>(num_mutations,0));

  // generate clone sequences based on clonal tree and mutations
  fprintf(stderr, "---\nNow generating clone genomes...\n");
  vector<Clone *> clones = tree.getVisibleNodes();
  //Clone root = *(tree.m_root);
  for (unsigned i=0; i<clones.size(); ++i)
    clones[i]->mutateGenome(ref_genome, mutations, variants, mutMatrix, clone2fn);

  // compile variants for output (only visible nodes in clone tree + root)
  //vector<Clone *> clones = tree.getVisibleNodes();
  vector<vector<short> > mat_mut_filt;
  vector<string> vec_labels;
  // add root
  vec_labels.push_back(tree.m_root->label);
  mat_mut_filt.push_back(mutMatrix[tree.m_root->index]);
  // add visible nodes
  for (unsigned i=0; i<clones.size(); i++) {
    vec_labels.push_back(clones[i]->label);
    mat_mut_filt.push_back(mutMatrix[clones[i]->index]);
  }

  // write clonal variants to file
  ofstream f_vcf;
  f_vcf.open("mutations.vcf");
  vario::writeVcf(ref_genome.records, variants, vec_labels, mat_mut_filt, f_vcf);
  f_vcf.close();

  // perform ADO
  if (ado_pct > 0)
    for (map<Clone*, string>::iterator it=clone2fn.begin(); it!=clone2fn.end(); ++it)
      seqio::simulateADO(it->second, ado_pct, ado_frag_len, random);

  return EXIT_SUCCESS;
}

/** Initialize random number generator
 */
boost::function<float()> initRandomNumberGenerator(long seed) {
  // init random number generator
  base_generator_type rng(seed);
  boost::uniform_real<> uni_dist(0, 1);
  boost::variate_generator<base_generator_type&, boost::uniform_real<> > uni(rng, uni_dist);

  boost::function<float()> f;
  f = boost::bind(uni_dist, rng);

  return f;
}

/** Parse command line arguments.
 * @return true: program can run normally, false: indication to stop
 */
bool parseArgs (int ac, char* av[], int& num_clones, vector<float>& freqs, int& num_mutations, int& num_transmuts, string& reference, string& ref_vcf, string& tree, bool verbose)
{
  std::stringstream ss;
  ss << std::endl << PROGRAM_NAME << std::endl << std::endl << "Available options:";

  namespace po = boost::program_options;

  // options only allowed on command line
  //po::options_description_generic();

  po::options_description desc(ss.str());
  desc.add_options()
    ("help,h", "print help message")
    ("clones,c", po::value<int>(&num_clones), "number of clones to simulate")
    ("freqs,f", po::value<vector<float> >(&freqs)->multitoken(), "clone relative frequencies")
    ("mutations,m", po::value<int>(&num_mutations), "total number of mutations")
    ("reference,r", po::value<string>(&reference), "reference sequence")
    ("reference-vcf,v", po::value<string>(&ref_vcf), "reference variants")
    ("initial-muts,i", po::value<int>(&num_transmuts)->default_value(1), "number of transforming mutations (separating healthy genome from first cancer genome)")
    ("tree,t", po::value<string>(&tree), "file containing user defined clone tree (Newick format)")
  ;

  po::variables_map var_map;

  try {
    po::store(po::parse_command_line(ac, av, desc), var_map);

    if (var_map.count("help") || ac == 1) {
      std::cerr << desc << std::endl;
      return false;
    }

    po::notify(var_map);  // might throw an error, so call after checking for "help"
  }
  catch (std::exception &e) {
    std::cerr << std::endl << e.what() << std::endl;
    std::cerr << desc << std::endl;
  }

  // check: num_clones > 0
  if (num_clones == 0 && !var_map.count("tree")) {
    fprintf(stderr, "\nArgumentError: Please specify clones>0 ('-c') or input tree ('-t'). Aborting...bye.\n");
    return false;
  }
  // check: #freqs == num_clones
  else if ((int)freqs.size() != num_clones) {
    fprintf(stderr, "\nArgumentError: Frequencies (%zu) need to match number of clones (%d).\n", freqs.size(), num_clones);
    return false;
  }
  // check: sum(freqs) == 1
  else {
    float sum_freqs = 0;
    for (unsigned int i=0; i<freqs.size(); i++) {
      sum_freqs += freqs[i];
    }
    if (sum_freqs != 1.0) {
      //fprintf(stderr, "\nArgumentError: Sum of frequencies (%.2f) needs to be 1.\n", sum_freqs);
      //return false;
    }
  }
  // check: num_mutations >= num_clones
  if (num_mutations < num_clones) {
    fprintf(stderr, "\nArgumentError: Number of mutations (%d) needs to be >= #clones (%d).\n", num_mutations, num_clones);
    return false;
  }
  // check: reference file exists
  if (var_map.count("reference")) {
    if (!fileExists(reference)) {
      fprintf(stderr, "\nArgumentError: File '%s' does not exist.\n", reference.c_str());
      return false;
    }
  }
  else {
    reference = "./data/chr7_114-115Mb.fa";
  }
  // check: reference VCF file exists
  if (var_map.count("reference-vcf")) {
    if (!fileExists(ref_vcf)) {
      fprintf(stderr, "\nArgumentError: File '%s' does not exist.\n", ref_vcf.c_str());
      return false;
    }
  }
  // check: was a clone tree provided by the user?
  if (var_map.count("tree")) {
    fprintf(stderr, "\nUser-defined clone tree was specified, parameter '-c/--clones' will be ignored\n");
    // check: does tree file exist?
    if (!fileExists(tree)) {
      fprintf(stderr, "\nArgumentError: File '%s' does not exist.\n", tree.c_str());
      return false;
    }
  }

  if (verbose) {
    fprintf(stderr, "---\n");
    fprintf(stderr, "Running with the following options:\n");
    if (var_map.count("tree")) {
      fprintf(stderr, "clone tree:\t%s\n", tree.c_str());
    } else {
      fprintf(stderr, "clones:\t\t%d\n", num_clones);
    }
    if (!var_map["freqs"].empty()) {
      fprintf(stderr, "frequencies:\t[ ");
      for (unsigned int i=0; i<freqs.size(); i++) {
        fprintf(stderr, "%.2f ", freqs[i]);
      }
      fprintf(stderr, "]\n");
    }
    fprintf(stderr, "mutations:\t%d\n", num_mutations);
    fprintf(stderr, "reference:\t%s\n", reference.c_str());
    if (var_map.count("reference-vcf")) {
      fprintf(stderr, "reference VCF:\t%s\n", ref_vcf.c_str());
    }
    fprintf(stderr, "---\n\n");
  }
  return true;
}

bool fileExists(string filename) {
  struct stat buffer;
  if (stat(filename.c_str(), &buffer)!=0) {
    return false;
  }
  return true;
}
