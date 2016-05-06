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
#include "yaml-cpp/yaml.h"

#define PROGRAM_NAME "CloniPhy 0.1"

using namespace std;
using seqio::SeqRecord;
using seqio::Genome;
using vario::Genotype;
using vario::Variant;
using evolution::SubstitutionModel;

boost::function<float()> initRandomNumberGenerator(long seed);
bool parseArgs (int ac, char* av[], YAML::Node&);
/*bool parseArgs (
  int ac, char* av[], string& conf, int& n_clones, std::vector<float>& freqs,
  int& n_mut, int& n_transmut, string& ref, string& ref_vcf, string& tree,
  bool verbose=true);*/
bool readConfig (string filename);
bool fileExists(string filename);

int main (int argc, char* argv[])
{
  // user params (defined in config file or command line)
  YAML::Node config = YAML::Node();
  bool args_ok = parseArgs(argc, argv, config);
  if (!args_ok) { return EXIT_FAILURE; }

  // params specified by user (in config file or command line)
  int num_clones = config["clones"].as<int>();
  int num_mutations = config["mutations"].as<int>();
  int num_transmuts = config["init-muts"].as<int>();
  vector<float> freqs = config["freqs"].as<vector<float> >();
  string reference = config["reference"].as<string>();
  string ref_vcf = config["reference-vcf"].as<string>();
  string tree_fn = config["tree"].as<string>();
  int verbosity = config["verbosity"].as<int>();
  long seed = config["seed"].as<long>();
  double titv = 1; // TODO: make this a user parameter (?)
  float ado_pct = 0.1; // TODO: make this a user parameter
  int ado_frag_len = 10000; // TODO: make this a user parameter

  // internal vars
  map<Clone*, string> clone2fn; // stores genome file names for each clone

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
      seqio::simulateADO(it->second, ref_genome.masked_length, ado_pct, ado_frag_len, random);

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
bool parseArgs (int ac, char* av[], YAML::Node& config)
  /*(int ac, char* av[], string& conf, int& num_clones, vector<float>& freqs,
  int& num_mutations, int& num_transmuts, string& reference, string& ref_vcf,
  string& tree, bool verbose)*/
{
  // default values
  int n_clones = 0;
  int n_mut = 0;
  int n_mut_init = 0;
  vector<float> freqs;
  string fn_ref = "";
  string fn_ref_vcf = "";
  string fn_tree = "";
  int verb = 1;
  long seed = time(0);

  stringstream ss;
  ss << endl << PROGRAM_NAME << endl << endl << "Available options";

  namespace po = boost::program_options;

  // options only allowed on command line
  //po::options_description_generic();

  po::options_description desc(ss.str());
  desc.add_options()
    ("help,h", "print help message")
    ("config,c", po::value<string>(), "config file")
    ("clones,n", po::value<int>(&n_clones), "number of clones to simulate")
    ("freqs,f", po::value<vector<float> >(&freqs)->multitoken(), "clone relative frequencies")
    ("mutations,m", po::value<int>(&n_mut), "total number of mutations")
    ("reference,r", po::value<string>(&fn_ref), "reference sequence")
    ("reference-vcf,v", po::value<string>(&fn_ref_vcf), "reference variants")
    ("init-muts,i", po::value<int>(&n_mut_init)->default_value(1), "number of transforming mutations (separating healthy genome from first cancer genome)")
    ("tree,t", po::value<string>(&fn_tree), "file containing user defined clone tree (Newick format)")
    ("verbosity,v", po::value<int>(&verb), "detail level of console output")
    ("seed,s", po::value<long>(&seed), "random seed")
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

  // check: config file exists
  if (var_map.count("config")) {
    string fn_config = var_map["config"].as<string>();
    if (!fileExists(fn_config)) {
      fprintf(stderr, "\nArgumentError: File '%s' does not exist.\n", fn_config.c_str());
      return false;
    }
    // initialize global configuration from config file
    config = YAML::LoadFile(fn_config);
  }

  // overwrite/set config params
  // (making sure parameters are set)
  if (var_map.count("clones") || !config["clones"]) {
    config["clones"] = n_clones;
  }
  n_clones = config["clones"].as<int>();
  if (var_map.count("freqs") || !config["freqs"]) {
    config["freqs"] = freqs;
  }
  freqs = config["freqs"].as<vector<float> >();
  if (var_map.count("mutations") || !config["mutations"]) {
    config["mutations"] = n_mut;
  }
  n_mut = config["mutations"].as<int>();
  if (var_map.count("init-muts") || !config["init-muts"]) {
    config["init-muts"] = n_mut_init;
  }
  n_mut_init = config["init-muts"].as<int>();
  if (var_map.count("reference") || !config["reference"]) {
    config["reference"] = fn_ref;
  }
  fn_ref = config["reference"].as<string>();
  if (var_map.count("reference-vcf") || !config["reference-vcf"]) {
    config["reference-vcf"] = fn_ref_vcf;
  }
  fn_ref = config["reference-vcf"].as<string>();
  if (var_map.count("tree") || !config["tree"]) {
    config["tree"] = fn_tree;
  }
  fn_tree = config["tree"].as<string>();
  if (var_map.count("verbosity") || !config["verbosity"]) {
    config["verbosity"] = verb;
  }
  verb = config["verbosity"].as<int>();
  if (var_map.count("seed") || !config["seed"]) {
    config["seed"] = seed;
  }
  seed = config["seed"].as<long>();


  // perform sanity checks

  // check: numClones > 0?
  if (n_clones == 0 && fn_tree.length() == 0) {
    fprintf(stderr, "\nArgumentError: Please specify clones ('-n') or input tree ('-t'). Aborting...bye.\n");
    return false;
  }
  // #freqs == num_clones?
  else if ((int)freqs.size() != n_clones) {
    fprintf(stderr, "\nArgumentError: Frequencies (%zu) need to match number of clones (%d).\n", freqs.size(), n_clones);
    return false;
  }
  // sum(freqs) == 1?
  else {
    float sum_freqs = 0;
    for (unsigned int i=0; i<freqs.size(); i++) {
      sum_freqs += freqs[i];
    }
    if (sum_freqs != 1.0) {
      fprintf(stderr, "\nArgumentError: Sum of frequencies (%.2f) needs to be 1.\n", sum_freqs);
      return false;
    }
  }
  // at least one mutation per clone?
  if (n_mut < n_clones) {
    fprintf(stderr, "\nArgumentError: Number of mutations (%d) needs to be >= #clones (%d).\n", n_mut, n_clones);
    return false;
  }
  // initial mutations do not exceed total mutations?
  if (n_mut_init > (n_mut-n_clones)) {
    fprintf(stderr, "\nArgumentError: Too many initial mutations (%d)\n -> can't be more than total mutations minus #clones (%d).\n", n_mut_init, n_mut-n_clones);
    return false;
  }
  // reference file exists?
  if (fn_ref.length()>0 && !fileExists(fn_ref)) {
    fprintf(stderr, "\nArgumentError: Reference genome file '%s' does not exist.\n", fn_ref.c_str());
    return false;
  }
  // reference VCF file exists?
  if (fn_ref_vcf.length()>0 && !fileExists(fn_ref_vcf)) {
    fprintf(stderr, "\nArgumentError: Reference VCF file '%s' does not exist.\n", fn_ref_vcf.c_str());
    return false;
  }
  // was a clone tree provided by the user?
  if (fn_tree.length()>0) {
    fprintf(stderr, "\nUser-defined clone tree was specified, parameter '-c/--clones' will be ignored\n");
    // check: does tree file exist?
    if (!fileExists(fn_tree)) {
      fprintf(stderr, "\nArgumentError: Tree file '%s' does not exist.\n", fn_tree.c_str());
      return false;
    }
  }

  if(config["verbosity"].as<int>() > 0) {
    fprintf(stderr, "---\n");
    fprintf(stderr, "Running with the following options:\n");
    if (fn_tree.length()>0) {
      fprintf(stderr, "clone tree:\t%s\n", fn_tree.c_str());
    } else {
      fprintf(stderr, "clones:\t\t%d\n", n_clones);
    }
    if (freqs.size()>0) {
      fprintf(stderr, "frequencies:\t[ ");
      for (unsigned int i=0; i<freqs.size(); i++) {
        fprintf(stderr, "%.2f ", freqs[i]);
      }
      fprintf(stderr, "]\n");
    }
    fprintf(stderr, "mutations:\t%d\n", n_mut);
    fprintf(stderr, "transforming mutations:\t%d\n", n_mut_init);
    fprintf(stderr, "reference:\t%s\n", fn_ref.c_str());
    if (fn_ref_vcf.length() > 0) {
      fprintf(stderr, "reference VCF:\t%s\n", fn_ref_vcf.c_str());
    }
    fprintf(stderr, "random seed:\t%ld\n", seed);
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
