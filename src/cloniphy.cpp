/**
 * Simulation of clonal DNA sequences based on the process of somatic evolution.
 *
 * @author: Harald Detering (harald.detering@gmail.com)
 *
 */

#include <boost/program_options.hpp>
namespace po = boost::program_options;
#include "simpleclonetree.hpp"
#include <ctime>

using namespace std;

bool parseArgs (int ac, char* av[], int *num_clones, std::vector<float> *freqs, int *num_mutations, bool verbose=true);

int main (int argc, char* argv[])
{
  // params specified by command line
  int num_clones;
  vector<float> freqs;
  int num_mutations;
  long int seed;

  fprintf(stderr, "This is CloniPhy. At the moment I am totally beta, so be patient :-)\n");
  bool args_ok = parseArgs(argc, argv, &num_clones, &freqs, &num_mutations);
  if (!args_ok) { return EXIT_FAILURE; }

  // specify random seed
  seed = 123456789;
  seed = time(NULL) + clock();

  SimpleCloneTree tree(num_clones, freqs);
  tree.generateRandomTopology(seed);

  fprintf(stderr, "\nNewick representation of generated tree:\n");
  Tree::printNewick(tree.getRoot(), cerr);
  fprintf(stderr, "\n");

  return EXIT_SUCCESS;
}

/** Parse command line arguments.
 * @return true if program is to be run normally, false to indicate to stop
 */
bool parseArgs (int ac, char* av[], int *num_clones, std::vector<float> *freqs, int *num_mutations, bool verbose)
{
  po::options_description desc("CanSim options");
  desc.add_options()
    ("help,h", "print help message")
    ("clones,c", po::value<int>(num_clones)->default_value(5), "number of clones to simulate")
    ("freqs,f", po::value<std::vector<float> >(freqs)->multitoken(), "clone relative frequencies")
    ("mutations,m", po::value<int>(num_mutations), "total number of mutations")
  ;

  po::variables_map var_map;
  po::store(po::parse_command_line(ac, av, desc), var_map);
  po::notify(var_map);

  if (var_map.count("help")) {
    std::cerr << desc << std::endl;
    return false;
  }
  // check: num_clones > 0
  if (*num_clones == 0) {
    fprintf(stderr, "\nArgumentError: 0 clones? Nothing to do here...bye.");
    return false;
  }
  // check: #freqs == num_clones
  else if ((int)freqs->size() != *num_clones) {
    fprintf(stderr, "\nArgumentError: Frequencies (%zu) need to match number of clones (%d).\n", freqs->size(), *num_clones);
    return false;
  }
  // check: sum(freqs) == 1
  else {
    float sum_freqs = 0;
    for (unsigned int i=0; i<freqs->size(); i++) {
      sum_freqs += (*freqs)[i];
    }
    if (sum_freqs != 1.0) {
      fprintf(stderr, "\nArgumentError: Sum of frequencies (%.2f) needs to be 1.\n", sum_freqs);
      return false;
    }
  }
  // check: num_mutations >= num_clones
  if (*num_mutations < *num_clones) {
    fprintf(stderr, "\nArgumentError: Number of mutations (%d) needs to be >= #clones (%d)\n", *num_mutations, *num_clones);
    return false;
  }

  if (verbose) {
    fprintf(stderr, "---\n");
    fprintf(stderr, "Running with the following options:\n");
    fprintf(stderr, "clones:\t\t%d\n", var_map["clones"].as<int>());
    if (!var_map["freqs"].empty()) {
      fprintf(stderr, "frequencies:\t[ ");
      for (unsigned int i=0; i<(*freqs).size(); i++) {
        fprintf(stderr, "%.2f ", (*freqs)[i]);
      }
      fprintf(stderr, "]\n");
    }
    fprintf(stderr, "mutations:\t%d\n", var_map["mutations"].as<int>());
    fprintf(stderr, "---\n\n");

  }
  return true;
}
