/**
 * Simulation of clonal DNA sequences based on the process of somatic evolution.
 *
 * @author: Harald Detering (harald.detering@gmail.com)
 *
 */
#include "seqio.hpp"
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
#include <set>
#include <sstream>
#include <string>
#include <sys/stat.h>

#define PROGRAM_NAME "CloniPhy 0.01"

using namespace std;

void generateMutations(std::vector<Mutation> &mutations, long ref_len, boost::function<float()>& random);
boost::function<float()> initRandomNumberGenerator(long seed);
bool parseArgs (int ac, char* av[], int& n_clones, std::vector<float>& freqs, int& n_mut, int& n_transmut, string& ref, bool verbose=true);

int main (int argc, char* argv[])
{
/*
std::string fn = "data/test.fa";
std::vector<SeqRecord> recs = SeqIO::readFasta(fn.c_str());
for (unsigned i=0; i<recs.size(); ++i) {
fprintf(stderr, ">%s\n%s...\n", recs[i].id.c_str(), recs[i].seq.substr(0,800).c_str());
}
return EXIT_SUCCESS;
*/
  // params specified by command line
  int num_clones = 0;
  int num_mutations = 0;
  int num_transmuts = 0;
  vector<float> freqs;
  string reference;
  long int seed;

  bool args_ok = parseArgs(argc, argv, num_clones, freqs, num_mutations, num_transmuts, reference);
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

  CoalescentCloneTree tree(num_clones, freqs);
  tree.generateRandomTopology(random);

  fprintf(stderr, "\nNewick representation of generated tree:\n");
  CoalescentCloneTree::printNewick(tree.getRoot(), cerr);
  fprintf(stderr, "\n");

  tree.evolve(num_mutations, num_transmuts, random);
  fprintf(stderr, "\nNewick representation of mutated tree:\n");
  CoalescentCloneTree::printNewick(tree.getRoot(), cerr);
  fprintf(stderr, "\n");

  fprintf(stderr, "Writing mutated tree to file 'clone_tree_mutated.dot'.\n");
  ofstream dotFileMut;
  dotFileMut.open("clone_tree_mutated.dot");
  tree.printDot(tree.getRoot(), dotFileMut);
  dotFileMut.close();

  tree.collapseZeroBranches(tree.getRoot());
  fprintf(stderr, "\nNewick representation of collapsed tree:\n");
  CoalescentCloneTree::printNewick(tree.getRoot(), cerr);
  fprintf(stderr, "\n");

  fprintf(stderr, "Writing collapsed tree to file 'clone_tree_collapsed.dot'.\n");
  ofstream dotFileCol;
  dotFileCol.open("clone_tree_collapsed.dot");
  tree.printDot(tree.getRoot(), dotFileCol);
  dotFileCol.close();

  // read reference sequence
  fprintf(stderr, "\nReading reference from file '%s'...", reference.c_str());
  std::vector<SeqRecord> refSeqs = SeqIO::readFasta(reference.c_str());
  unsigned i;
  unsigned refLen = 0;
  std::vector<long> cumStart; // cumulative start position of each sequence;
  for (i=0; i<refSeqs.size(); ++i) {
    cumStart.push_back(refLen);
    refLen += refSeqs[i].seq.size();
  }
  fprintf(stderr, "read (%u bp in %u sequences).\n", refLen, i);
  fprintf(stderr, "generating FASTA index.\n");
  SeqIO::indexFasta(reference.c_str());

  // generate relative positions of mutations;
  std::vector<Mutation> mutations(num_mutations);
  generateMutations(mutations, refLen, random);

  fprintf(stderr, "\nTotal set of mutations (bp, offset):\n");
  for (int i=0; i<num_mutations; i++) {
    fprintf(stderr, "%ld\t%d\n", mutations[i].absPos, mutations[i].offset);
  }
  // generate clone sequences based on clonal tree and mutations
  std::vector<Clone *> leafs = tree.getLeafs();
  for (std::vector<Clone *>::iterator i=leafs.begin(); i!=leafs.end(); ++i) {
    Clone c = **i;
//fprintf(stderr, "<Clone(label=%d)>, %u mutations\n", c.label, c.m_vecMutations.size());
    std::vector<SeqRecord> cloneGenome = refSeqs; // TODO: check memory footprint
    for (unsigned i=0; i<cloneGenome.size(); ++i) {
      std::string seqId = boost::str(boost::format("clone%02d_seq%d") % c.label % (i+1));
      cloneGenome[i].id = seqId;
    }
    c.mutateGenome(cloneGenome, cumStart, mutations);
    std::string filename = boost::str(boost::format("clone%02d.fa") % c.label);
    ofstream outfile;
    outfile.open(filename.c_str());
    SeqIO::writeFasta(cloneGenome, outfile);
  }

  return EXIT_SUCCESS;
}

/** Generate random mutations out of thin air. */
void generateMutations(std::vector<Mutation> &mutations, long ref_len, boost::function<float()>& random) {
  std::set<long> mutPositions; // remember mutated positions (enforce infinite sites model)
  for (std::vector<Mutation>::iterator m=mutations.begin(); m!=mutations.end(); ++m) {
    float rel_pos = random();
    long abs_pos = rel_pos * ref_len;
    // enforce infinite sites (no position can be mutated more than once)
    while (mutPositions.count(abs_pos)!=0) {
      abs_pos = (abs_pos+1) % ref_len; // wrap around at end of reference
    }
    int offset = (random()*3)+1;
fprintf(stderr, "<Mutation(absPos=%ld,offset=%d)>\n", abs_pos, offset);
    m->absPos = abs_pos;
    m->offset = offset;
  }
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
bool parseArgs (int ac, char* av[], int& num_clones, std::vector<float>& freqs, int& num_mutations, int& num_transmuts, std::string& reference, bool verbose)
{
  std::stringstream ss;
  ss << std::endl << PROGRAM_NAME << std::endl << std::endl << "Available options:";

  namespace po = boost::program_options;
  po::options_description desc(ss.str());
  desc.add_options()
    ("help,h", "print help message")
    ("clones,c", po::value<int>(&num_clones), "number of clones to simulate")
    ("freqs,f", po::value<std::vector<float> >(&freqs)->multitoken(), "clone relative frequencies")
    ("mutations,m", po::value<int>(&num_mutations), "total number of mutations")
    ("reference,r", po::value<std::string>(&reference), "reference sequence")
    ("trans-muts,t", po::value<int>(&num_transmuts)->default_value(1), "number of transforming mutations (separating healthy genome from first cancer genome)")
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
  if (num_clones == 0) {
    fprintf(stderr, "\nArgumentError: 0 clones? Nothing to do here...bye.\n");
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
      fprintf(stderr, "\nArgumentError: Sum of frequencies (%.2f) needs to be 1.\n", sum_freqs);
      return false;
    }
  }
  // check: num_mutations >= num_clones
  if (num_mutations < num_clones) {
    fprintf(stderr, "\nArgumentError: Number of mutations (%d) needs to be >= #clones (%d).\n", num_mutations, num_clones);
    return false;
  }
  // check: reference file exists
  if (var_map.count("reference")) {
    struct stat buffer;
    if (stat(reference.c_str(), &buffer)!=0) {
      fprintf(stderr, "\nArgumentError: File '%s' does not exist.\n", reference.c_str());
      return false;
    }
  }
  else {
    reference = "./data/chr7_114-115Mb.fa";
  }

  if (verbose) {
    fprintf(stderr, "---\n");
    fprintf(stderr, "Running with the following options:\n");
    fprintf(stderr, "clones:\t\t%d\n", var_map["clones"].as<int>());
    if (!var_map["freqs"].empty()) {
      fprintf(stderr, "frequencies:\t[ ");
      for (unsigned int i=0; i<freqs.size(); i++) {
        fprintf(stderr, "%.2f ", freqs[i]);
      }
      fprintf(stderr, "]\n");
    }
    fprintf(stderr, "mutations:\t%d\n", var_map["mutations"].as<int>());
    fprintf(stderr, "---\n\n");
  }
  return true;
}
