#include "ConfigStore.hpp"

using namespace std;

namespace config {

SampleConfig::SampleConfig(YAML::Node row) {
  this->m_label = row[0].as<string>();
  this->m_vec_prevalence.clear();
  for (int i=1; i<row.size(); i++) {
    this->m_vec_prevalence.push_back(row[i].as<double>());
  }
}

ConfigStore::ConfigStore() {
  _config = YAML::Node();
}

/** Parse command line arguments.
 * @return true: program can run normally, false: indication to stop
 */
bool ConfigStore::parseArgs (int ac, char* av[])
{
  // default values
  int n_clones = 0;
  int n_mut = 0;
  int n_mut_init = 0;
  string model = "JC";
  string fn_bam_input = "";
  string fn_mut_ref_vcf = "";
  string fn_mut_som_vcf = "";
  string fn_mut_sig = "resources/signatures_probabilities.txt";
  string fn_ref_fa = "";
  string fn_tree = "";
  string pfx_out = "clonesim";
  int verb = 1;
  long seed = time(NULL) + clock();

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
    ("mutations,m", po::value<int>(&n_mut), "total number of mutations")
    ("reference,r", po::value<string>(&fn_ref_fa), "reference sequence")
    ("reference-vcf,v", po::value<string>(&fn_mut_ref_vcf), "reference variants")
    ("init-muts,i", po::value<int>(&n_mut_init), "number of transforming mutations (separating healthy genome from first cancer genome)")
    ("tree,t", po::value<string>(&fn_tree), "file containing user defined clone tree (Newick format)")
    ("out,o", po::value<string>(&pfx_out), "prefix for output files")
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
  if (var_map.count("mutations") || !_config["mutations"]) {
    _config["mutations"] = n_mut;
  }
  n_mut = _config["mutations"].as<int>();
  if (var_map.count("init-muts") || !_config["init-muts"]) {
    _config["init-muts"] = n_mut_init;
  }
  n_mut_init = _config["init-muts"].as<int>();
  if (var_map.count("reference") || !_config["reference"]) {
    _config["reference"] = fn_ref_fa;
  }
  if (_config["mut-sig-file"]) {
    fn_mut_sig = _config["mut-sig-file"].as<string>();
  } else {
    _config["mut-sig-file"] = fn_mut_sig;
  }
  if (_config["bam-input"]) {
    fn_bam_input = _config["bam-input"].as<string>();
  } else {
    _config["bam-input"] = fn_bam_input;
  }
  fn_ref_fa = _config["reference"].as<string>();
  if (var_map.count("reference-vcf") || !_config["reference-vcf"]) {
    _config["reference-vcf"] = fn_mut_ref_vcf;
  }
  fn_mut_ref_vcf = _config["reference-vcf"].as<string>();
  if (var_map.count("mut-som-vcf") || !_config["mut-som-vcf"]) {
    _config["mut-som-vcf"] = fn_mut_som_vcf;
  }
  fn_mut_som_vcf = _config["mut-som-vcf"].as<string>();
  if (var_map.count("tree") || !_config["tree"]) {
    _config["tree"] = fn_tree;
  }
  fn_tree = _config["tree"].as<string>();
  if (var_map.count("out") || !_config["out"]) {
    _config["out"] = pfx_out;
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

  // at least one mutation per clone?
  if (n_mut < n_clones) {
    if (n_mut > 0) {
      fprintf(stderr, "\nArgumentError: Number of mutations (%d) needs to be >= #clones (%d).\n", n_mut, n_clones);
      return false;
    } else if (n_mut == 0 && fn_tree.length()==0) {
      fprintf(stderr, "\nArgumentError: When setting number of mutations = 0, a tree file needs to be specified. (Branch lengths -> #mutations)\n");
      return false;
    }
  }
  // initial mutations do not exceed total mutations?
  //if (n_mut_init > (n_mut-n_clones)) {
  //  fprintf(stderr, "\nArgumentError: Too many initial mutations (%d)\n -> can't be more than total mutations minus #clones (%d).\n", n_mut_init, n_mut-n_clones);
  //  return false;
  //}
  // input BAM file exists?
  if (fn_bam_input.length()>0 && !fileExists(fn_bam_input)) {
    fprintf(stderr, "\nArgumentError: Input BAM file '%s' does not exist.\n", fn_bam_input.c_str());
    return false;
  }
  // reference file exists?
  if (fn_ref_fa.length()>0 && !fileExists(fn_ref_fa)) {
    fprintf(stderr, "\nArgumentError: Reference genome file '%s' does not exist.\n", fn_ref_fa.c_str());
    return false;
  }
  // reference VCF file exists?
  if (fn_mut_ref_vcf.length()>0 && !fileExists(fn_mut_ref_vcf)) {
    fprintf(stderr, "\nArgumentError: Reference VCF file '%s' does not exist.\n", fn_mut_ref_vcf.c_str());
    return false;
  }
  // somatic VCF file exists?
  if (fn_mut_som_vcf.length()>0 && !fileExists(fn_mut_som_vcf)) {
    fprintf(stderr, "\nArgumentError: Somatic VCF file '%s' does not exist.\n", fn_mut_som_vcf.c_str());
    return false;
  }
  // was a clone tree provided by the user?
  if (fn_tree.length()>0) {
    //fprintf(stderr, "\nUser-defined clone tree was specified, parameter '-c/--clones' will be ignored\n");
    // check: does tree file exist?
    if (!fileExists(fn_tree)) {
      fprintf(stderr, "\nArgumentError: Tree file '%s' does not exist.\n", fn_tree.c_str());
      return false;
    }
  }

  // check germline evolutionary model params
  if (!_config["model"]) {
    //fprintf(stderr, "\n[INFO] Missing evolutionary model - assuming '%s'.\n", model.c_str());
    YAML::Node node = YAML::Load(model);
    _config["model"] = node;
  } else {
    model = _config["model"].as<string>();
  }
  if (model != "JC" && !_config["model-params"]) {
    fprintf(stderr, "\nArgumentError: Evolutionary models other than 'JC' require parameters (set param 'model-params' in config file.)\n");
    return false;
  }
  if (model == "F81" || model == "HKY") {
    if (!_config["model-params"]["nucFreq"]) {
      fprintf(stderr, "\nArgumentError: Model '%s' requires nucleotide frequencies (set param 'model-params':'nucFreq' in config file.)\n", model.c_str());
      return false;
    } else if (_config["model-params"]["nucFreq"].size() != 4) {
      fprintf(stderr, "\nArgumentError: Parameter 'model-params':'nucFreq' must contain exactly 4 values.\n");
      return false;
    }
  }
  if (model == "K80" || model == "HKY") {
    if (!_config["model-params"]["kappa"]) {
      fprintf(stderr, "\nArgumentError: Model '%s' requires transition-transversion ratio (set param 'model-params':'kappa' in config file.)\n", model.c_str());
      return false;
    }
  }
  if (model == "matrix") {
    vector<string> names = { "Qa", "Qc", "Qg", "Qt" };
    for (string n : names) {
      if (!_config["model-params"][n] ) {
        fprintf(stderr, "\nArgumentError: Model 'matrix' requires substitution frequencies (set param 'model-params':'Qa','Qc','Qg','Qt' in config file.)");
        fprintf(stderr, "\n               Missing parameter: '%s'.)\n", n.c_str());
        return false;
      } else if (_config["model-params"][n].as<vector<double>>().size() != 4 ) {
        fprintf(stderr, "\nArgumentError: Substitution frequencies vectors must contain exactly 4 real values.");
        fprintf(stderr, "\n               Violating parameter: 'model-params':'%s'.)\n", n.c_str());
        return false;
      }
    }
  }

  // check somatic evolution parameters
  if (!fileExists(fn_mut_sig)) {
    fprintf(stderr, "\nArgumentError: Mutation signature file '%s' does not exist.\n", fn_mut_sig.c_str());
    return false;
  }
  if (!_config["mut-signatures"]) {
    fprintf(stderr, "\n[INFO] Missing somatic mutation signatures - assuming signature 1.\n");
    YAML::Node node = YAML::Load("1: 1.0");
    _config["mut-signatures"] = node;
  }

  // does sampling matrix have the expected number of rows (clones + 1)?
  if (!_config["samples"] || _config["samples"].size() == 0) {
    fprintf(stderr, "\nArgumentError: Missing sampling matrix - 'samples' param in config file needed.\n");
    return false;
  } else if (fn_tree.empty()) {
    double row_sum;
    for (auto row_sample : _config["samples"]) {
      if (row_sample.size() != (n_clones + 1)) {
        fprintf(stderr, "\nArgumentError: Columns in sampling matrix must be clones+1 (violated in '%s').\n", row_sample[0].as<string>().c_str());
        return false;
      }
      // check: row sum <= 1?
      row_sum = 0.0;
      for (size_t i=1; i<row_sample.size(); i++)
        row_sum += row_sample[i].as<double>();
      if (row_sum > 1) {
        fprintf(stderr, "\nArgumentError: Row sums in sampling matrix must <=1 (violated in '%s').\n", row_sample[0].as<string>().c_str());
        return false;
      }
      this->m_vec_samples.push_back(SampleConfig(row_sample));
    }
  }

  if(_config["verbosity"].as<int>() > 0) {
    fprintf(stderr, "---\n");
    fprintf(stderr, "Running with the following options:\n");
    fprintf(stderr, "==================================\n");
    fprintf(stderr, "  random seed:\t\t%ld\n", seed);
    if (fn_tree.length()>0) {
      fprintf(stderr, "  clone tree:\t\t%s\n", fn_tree.c_str());
    } else {
      fprintf(stderr, "  clones:\t\t%d\n", n_clones);
    }
    fprintf(stderr, "----------------------------------\n");
    fprintf(stderr, "Reference genome:\n");
    fprintf(stderr, "----------------------------------\n");
    if (fn_ref_fa.length()>0) {
      fprintf(stderr, "  input file:\t\t%s\n", fn_ref_fa.c_str());
    } else {
      fprintf(stderr, "  generate in-silico:\tyes\n");
      fprintf(stderr, "  ref seqs:\t\t%d\n", this->getValue<int>("ref-seq-num"));
      fprintf(stderr, "  seq len:\t\t%d (+/-%d)\n", this->getValue<int>("ref-seq-len-mean"), this->getValue<int>("ref-seq-len-mean"));
    }
    fprintf(stderr, "----------------------------------\n");
    fprintf(stderr, "Germline mutations:\n");
    fprintf(stderr, "----------------------------------\n");
    if (fn_mut_ref_vcf.length() > 0) {
      fprintf(stderr, "reference VCF:\t%s\n", fn_mut_ref_vcf.c_str());
    } else {
      // TODO: print germline mutation params
    }
    fprintf(stderr, "----------------------------------\n");
    fprintf(stderr, "Sequencing data\n");
    fprintf(stderr, "----------------------------------\n");
    if (fn_bam_input.length()>0) {
      fprintf(stderr, "  simulate reads:\tno\n");
      fprintf(stderr, "  input reads:\t%s\n", fn_bam_input.c_str());
    } else {
      fprintf(stderr, "  simulate reads:\tyes\n");
      fprintf(stderr, "  ref coverage:\t\t%d\n", this->getValue<int>("seq-coverage"));
      fprintf(stderr, "  seq read length:\t%d\n", this->getValue<int>("seq-read-len"));
      fprintf(stderr, "  seq insert size:\t%d (+-%d)\n", this->getValue<int>("seq-frag-len-mean"), this->getValue<int>("seq-frag-len-sd"));
      fprintf(stderr, "  simulator:\t\t%s\n", this->getValue<string>("seq-art-path").c_str());
    }
    fprintf(stderr, "----------------------------------\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "----------------------------------\n");
    fprintf(stderr, "mutations:\t%d\n", n_mut);
    fprintf(stderr, "transforming mutations:\t%d\n", n_mut_init);
    fprintf(stderr, "evolutionary model:\t%s\n", model.c_str());
    if (_config["samples"]) {
      fprintf(stderr, "Sampling scheme:\n");
      map<string, vector<double>> sample_mtx = this->getMatrix<double>("samples");
      string s_sampling = stringio::printMatrix(sample_mtx);
      fprintf(stderr, "%s", s_sampling.c_str());
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

} /* namespace config */
