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
#include <sstream>
#include <string>
#include <sys/stat.h>

#define PROGRAM_NAME "CloniPhy 0.1"

using namespace std;

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
  }
  else {
    //CoalescentCloneTree tree(num_clones, freqs);
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
  vector<SeqRecord> ref_seqs = SeqIO::readFasta(reference.c_str());
  unsigned i;
  unsigned long ref_len = 0;
  vector<unsigned long> cumStart; // cumulative start position of each sequence;
  for (i=0; i<ref_seqs.size(); ++i) {
    cumStart.push_back(ref_len);
    ref_len += ref_seqs[i].seq.size();
  }
  fprintf(stderr, "read (%lu bp in %u sequences).\n", ref_len, i);
  if (ref_vcf.size() > 0) {
    fprintf(stderr, "applying reference variants (from %s).\n", ref_vcf.c_str());
    vector<Variant> ref_variants;
    vector<vector<Genotype > > ref_gt_matrix;
    VarIO::readVcf(ref_vcf, ref_variants, ref_gt_matrix);
    VarIO::applyVariants(ref_seqs, ref_variants, ref_gt_matrix[0]);
  }
  else {
    // duplicate genome (all loci homozygous reference)
    unsigned num_seqs = ref_seqs.size();
    for (unsigned i=0; i<num_seqs; ++i) {
      SeqRecord orig = ref_seqs[i];
      SeqRecord *dupl = new SeqRecord(orig.id, orig.description, orig.seq);
      ref_seqs.push_back(*dupl);
      ref_seqs[i].id += "_m";
      ref_seqs[num_seqs+i].id += "_p";
    }
  }
  // write "healthy" genome to file
  ofstream f_fasta;
  f_fasta.open("healthy_genome.fa");
  SeqIO::writeFasta(ref_seqs, f_fasta);
  f_fasta.close();

  //fprintf(stderr, "generating FASTA index.\n");
  //SeqIO::indexFasta(reference.c_str());

  // generate mutations (absolute position + nucleotide shift + chr copy)
  vector<Mutation> mutations = VarIO::generateMutations(num_mutations, ref_len, random);

  fprintf(stderr, "\nTotal set of mutations (bp, offset):\n");
  for (int i=0; i<num_mutations; i++) {
    fprintf(stderr, "%ld\t%d\n", mutations[i].absPos, mutations[i].offset);
  }

  // initialize mutation matrix
  vector<std::vector<short> > mutMatrix(num_clones+1, std::vector<short>(num_mutations,0));

  // generate clone sequences based on clonal tree and mutations
  vector<Clone *> clones = tree.getVisibleNodes();
  for (unsigned i=0; i<clones.size(); ++i) {
    Clone c = *(clones[i]);
cerr << c << ", " << c.m_vecMutations.size() << " mutations" << endl;
    vector<SeqRecord> clone_genome = ref_seqs; // TODO: check memory footprint
    // assign clone id (either by label or index)
    string clone_id = "";
    if (c.label.size()>0) {
      clone_id = c.label;
    } else {
      clone_id = boost::str(boost::format("clone_%02d") % c.index);
    }
    // rename sequences including clone id
    for (unsigned j=0; j<clone_genome.size(); ++j) {
      clone_genome[j].id += '-' + clone_id;
    }
    c.mutateGenome(clone_genome, cumStart, mutations, mutMatrix[i]);
    string filename = boost::str(boost::format("%s_genome.fa") % clone_id);
    ofstream outfile;
    outfile.open(filename.c_str());
    SeqIO::writeFasta(clone_genome, outfile);
  }

  ofstream f_vcf;
  f_vcf.open("mutations.vcf");
  VarIO::writeVcf(ref_seqs, mutations, mutMatrix, f_vcf);
  f_vcf.close();

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
