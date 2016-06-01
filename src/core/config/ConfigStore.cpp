#include "ConfigStore.hpp"

using namespace std;

namespace config {

ConfigStore::ConfigStore() {
  _config = YAML::Node();
}

/** Parse command line arguments.
 * @return true: program can run normally, false: indication to stop
 */
bool ConfigStore::parseArgs (int ac, char* av[])
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
    _config = YAML::LoadFile(fn_config);
  }

  // overwrite/set config params
  // (making sure parameters are set)
  if (var_map.count("clones") || !_config["clones"]) {
    _config["clones"] = n_clones;
  }
  n_clones = _config["clones"].as<int>();
  if (var_map.count("freqs") || !_config["freqs"]) {
    _config["freqs"] = freqs;
  }
  freqs = _config["freqs"].as<vector<float> >();
  if (var_map.count("mutations") || !_config["mutations"]) {
    _config["mutations"] = n_mut;
  }
  n_mut = _config["mutations"].as<int>();
  if (var_map.count("init-muts") || !_config["init-muts"]) {
    _config["init-muts"] = n_mut_init;
  }
  n_mut_init = _config["init-muts"].as<int>();
  if (var_map.count("reference") || !_config["reference"]) {
    _config["reference"] = fn_ref;
  }
  fn_ref = _config["reference"].as<string>();
  if (var_map.count("reference-vcf") || !_config["reference-vcf"]) {
    _config["reference-vcf"] = fn_ref_vcf;
  }
  fn_ref = _config["reference-vcf"].as<string>();
  if (var_map.count("tree") || !_config["tree"]) {
    _config["tree"] = fn_tree;
  }
  fn_tree = _config["tree"].as<string>();
  if (var_map.count("verbosity") || !_config["verbosity"]) {
    _config["verbosity"] = verb;
  }
  verb = _config["verbosity"].as<int>();
  if (var_map.count("seed") || !_config["seed"]) {
    _config["seed"] = seed;
  }
  seed = _config["seed"].as<long>();


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

  if(_config["verbosity"].as<int>() > 0) {
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

}
