#include "ConfigStore.hpp"
#include <gitversion/version.h>

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
  int n_mut_somatic = 0;
  int n_mut_trunk = 0;
  int n_mut_gl = 0;
  double mut_som_cnv_ratio = 0.0;
  string mut_gl_model = "JC";
  string dir_out = "output";
  string fn_bam_input = "";
  string fn_mut_gl_vcf = "";
  string fn_mut_som_vcf = "";
  string fn_mut_som_sig = "resources/signatures_probabilities.txt";
  string fn_ref_fa = "";
  string fn_tree = "";
  bool seq_read_gen = false;
  int seq_rc_min = 1;
  bool seq_use_vaf = false;
  bool do_reuse_reads = false;
  bool do_fq_out = true;
  bool do_sam_out = true;
  int verb = 1;
  long seed = time(NULL) + clock();

  stringstream ss;
  ss << endl << PROGRAM_NAME << " " << version::GIT_TAG_NAME << endl << endl;
  ss << "Available options";

  namespace po = boost::program_options;

  // options only allowed on command line
  //po::options_description_generic();

  po::options_description desc(ss.str());
  desc.add_options()
    ("version,v", "print version string")
    ("help,h", "print help message")
    ("config,c", po::value<string>(), "config file")
    ("clones,n", po::value<int>(&n_clones), "number of clones to simulate")
    ("mut-som-num,m", po::value<int>(&n_mut_somatic), "total number of mutations")
    ("reference,r", po::value<string>(&fn_ref_fa), "reference sequence")
    ("mut-gl-vcf,v", po::value<string>(&fn_mut_gl_vcf), "germline variants")
    ("mut-som-trunk,i", po::value<int>(&n_mut_trunk), "number of transforming mutations (separating healthy genome from first cancer genome)")
    ("tree,t", po::value<string>(&fn_tree), "file containing user defined clone tree (Newick format)")
    ("out-dir,o", po::value<string>(&dir_out), "output directory")
    ("verbosity,v", po::value<int>(&verb), "detail level of console output")
    ("seed,s", po::value<long>(&seed), "random seed")
  ;

  po::variables_map var_map;

  try {
    po::store(po::parse_command_line(ac, av, desc), var_map);

    if (var_map.count("version")) {
      std::cerr << PROGRAM_NAME << " " << version::GIT_TAG_NAME << endl;
      return false;
    }

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

  // number of clones to simulate
  if (var_map.count("clones") || !_config["clones"]) {
    _config["clones"] = n_clones;
  }
  n_clones = _config["clones"].as<int>();
  // number of somatic mutations
  if (var_map.count("mut-som-num") || !_config["mut-som-num"]) {
    _config["mut-som-num"] = n_mut_somatic;
  }
  n_mut_somatic = _config["mut-som-num"].as<int>();
  // number of somatic mutations assigned to trunk of clone tree
  if (var_map.count("mut-som-trunk") || !_config["mut-som-trunk"]) {
    _config["mut-som-trunk"] = n_mut_trunk;
  }
  n_mut_trunk = _config["mut-som-trunk"].as<int>();
  // path to reference FASTA
  if (var_map.count("reference") || !_config["reference"]) {
    _config["reference"] = fn_ref_fa;
  }
  fn_ref_fa = _config["reference"].as<string>();
  // path to somatic mutation signature file
  if (_config["mut-som-sig-file"]) {
    fn_mut_som_sig = _config["mut-som-sig-file"].as<string>();
  } else {
    _config["mut-som-sig-file"] = fn_mut_som_sig;
  }
  if (_config["bam-input"]) {
    fn_bam_input = _config["bam-input"].as<string>();
  } else {
    _config["bam-input"] = fn_bam_input;
  }
  // number of germline mutations
  if (!_config["mut-gl-num"]) {
    _config["mut-gl-num"] = n_mut_gl;
  }
  n_mut_gl = _config["mut-gl-num"].as<int>();
  // path to germline variants input VCF
  if (var_map.count("mut-gl-vcf") || !_config["mut-gl-vcf"]) {
    _config["mut-gl-vcf"] = fn_mut_gl_vcf;
  }
  fn_mut_gl_vcf = _config["mut-gl-vcf"].as<string>();
  // path to somatic variants input VCF
  if (var_map.count("mut-som-vcf") || !_config["mut-som-vcf"]) {
    _config["mut-som-vcf"] = fn_mut_som_vcf;
  }
  fn_mut_som_vcf = _config["mut-som-vcf"].as<string>();
  // fraction of CNV events among all somatic mutations
  if (!_config["mut-som-cnv-ratio"]) {
    _config["mut-som-cnv-ratio"] = mut_som_cnv_ratio;
  }
  mut_som_cnv_ratio = _config["mut-som-cnv-ratio"].as<double>();
  // path to clone tree input file (Newick)
  if (var_map.count("tree") || !_config["tree"]) {
    _config["tree"] = fn_tree;
  }
  fn_tree = _config["tree"].as<string>();
  // prefix to use for output files
  if (var_map.count("out-dir") || !_config["out-dir"]) {
    _config["out-dir"] = dir_out;
  }
  dir_out = _config["out-dir"].as<string>();
  // how chatty should status messages be?
  if (var_map.count("verbosity") || !_config["verbosity"]) {
    _config["verbosity"] = verb;
  }
  verb = _config["verbosity"].as<int>();
  // seed for random number generator
  if (var_map.count("seed") || !_config["seed"]) {
    _config["seed"] = seed;
  }
  seed = _config["seed"].as<long>();

  //---------------------------------------------------------------------------
  // sequencing-related params
  //---------------------------------------------------------------------------

  // whether to generate sequencing reads (if not, generate readcounts directly)
  if (!_config["seq-read-gen"]) {
    _config["seq-read-gen"] = seq_read_gen;
  }
  seq_read_gen = _config["seq-read-gen"].as<bool>();
  // whether to use VAFs to spike in mutations (if not, assign read pairs to segment copies)
  if (!_config["seq-use-vaf"]) {
    _config["seq-use-vaf"] = seq_use_vaf;
  }
  seq_use_vaf = _config["seq-use-vaf"].as<bool>();
  // bit to indicate if reads of normal sample are be reused for tumor samples
  if (!_config["seq-reuse-reads"]) {
      _config["seq-reuse-reads"] = do_reuse_reads;
  }
  do_reuse_reads = _config["seq-reuse-reads"].as<bool>();
  // bit to indicate if sequence simulator should keep FASTQ files
  if (!_config["seq-fq-out"]) {
      _config["seq-fq-out"] = do_fq_out;
  }
  do_fq_out = _config["seq-fq-out"].as<bool>();
  // bit to indicate if sequence simulator should generate SAM files
  if (!_config["seq-sam-out"]) {
      _config["seq-sam-out"] = do_sam_out;
  }
  do_fq_out = _config["seq-sam-out"].as<bool>();
  // when generating read counts, minimum ALT read count for which to output a VCF line
  if (!_config["seq-rc-min"]) {
    _config["seq-rc-min"] = seq_rc_min;
  }
  seq_rc_min = _config["seq-rc-min"].as<int>();
  
  //---------------------------------------------------------------------------
  // perform sanity checks
  //---------------------------------------------------------------------------

  // at least one mutation per clone?
  if (n_mut_somatic < n_clones) {
    if (n_mut_somatic > 0) {
      fprintf(stderr, "\nArgumentError: Number of mutations (%d) needs to be >= #clones (%d).\n", n_mut_somatic, n_clones);
      return false;
    } else if (n_mut_somatic == 0 && fn_tree.length()==0) {
      fprintf(stderr, "\nArgumentError: When setting number of mutations = 0, a tree file needs to be specified. (Branch lengths -> #mutations)\n");
      return false;
    }
  }
  // initial mutations do not exceed total mutations?
  //if (n_mut_trunk > (n_mut-n_clones)) {
  //  fprintf(stderr, "\nArgumentError: Too many initial mutations (%d)\n -> can't be more than total mutations minus #clones (%d).\n", n_mut_trunk, n_mut-n_clones);
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
  if (fn_mut_gl_vcf.length()>0 && !fileExists(fn_mut_gl_vcf)) {
    fprintf(stderr, "\nArgumentError: Germline VCF file '%s' does not exist.\n", fn_mut_gl_vcf.c_str());
    return false;
  }
  // somatic VCF file exists?
  if (fn_mut_som_vcf.length()>0 && !fileExists(fn_mut_som_vcf)) {
    fprintf(stderr, "\nArgumentError: Somatic VCF file '%s' does not exist.\n", fn_mut_som_vcf.c_str());
    return false;
  }
  // somatic mutation signature file exists?
  if (fn_mut_som_sig.length()>0 && !fileExists(fn_mut_som_sig)) {
    fprintf(stderr, "\nArgumentError: Somatic mutation signature file '%s' does not exist.\n", fn_mut_som_sig.c_str());
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
  if (!_config["mut-gl-model"]) {
    //fprintf(stderr, "\n[INFO] Missing evolutionary model - assuming '%s'.\n", model.c_str());
    YAML::Node node = YAML::Load(mut_gl_model);
    _config["mut-gl-model"] = node;
  } else {
    mut_gl_model = _config["mut-gl-model"].as<string>();
  }
  if (mut_gl_model != "JC" && !_config["mut-gl-model-params"]) {
    fprintf(stderr, "\nArgumentError: Evolutionary models other than 'JC' require parameters (set param 'mut-gl-model-params' in config file.)\n");
    return false;
  }
  if (mut_gl_model == "F81" || mut_gl_model == "HKY") {
    if (!_config["mut-gl-model-params"]["nucFreq"]) {
      fprintf(stderr, "\nArgumentError: Model '%s' requires nucleotide frequencies (set param 'mut-gl-model-params':'nucFreq' in config file.)\n", mut_gl_model.c_str());
      return false;
    } else if (_config["mut-gl-model-params"]["nucFreq"].size() != 4) {
      fprintf(stderr, "\nArgumentError: Parameter 'mut-gl-model-params':'nucFreq' must contain exactly 4 values.\n");
      return false;
    }
  }
  if (mut_gl_model == "K80" || mut_gl_model == "HKY") {
    if (!_config["mut-gl-model-params"]["kappa"]) {
      fprintf(stderr, "\nArgumentError: Model '%s' requires transition-transversion ratio (set param 'mut-gl-model-params':'kappa' in config file.)\n", mut_gl_model.c_str());
      return false;
    }
  }
  if (mut_gl_model == "matrix") {
    vector<string> names = { "Qa", "Qc", "Qg", "Qt" };
    for (string n : names) {
      if (!_config["mut-gl-model-params"][n] ) {
        fprintf(stderr, "\nArgumentError: Model 'matrix' requires substitution frequencies (set param 'mut-gl-model-params':'Qa','Qc','Qg','Qt' in config file.)");
        fprintf(stderr, "\n               Missing parameter: '%s'.)\n", n.c_str());
        return false;
      } else if (_config["mut-gl-model-params"][n].as<vector<double>>().size() != 4 ) {
        fprintf(stderr, "\nArgumentError: Substitution frequencies vectors must contain exactly 4 real values.");
        fprintf(stderr, "\n               Violating parameter: 'mut-gl-model-params':'%s'.)\n", n.c_str());
        return false;
      }
    }
  }

  // check somatic evolution parameters
  if (!fileExists(fn_mut_som_sig)) {
    fprintf(stderr, "\nArgumentError: Mutation signature file '%s' does not exist.\n", fn_mut_som_sig.c_str());
    return false;
  }
  if (!_config["mut-som-sig-mix"]) {
    fprintf(stderr, "\n[INFO] Missing somatic mutation signatures - assuming signature 1.\n");
    YAML::Node node = YAML::Load("{'Signature 1': 1.0}");
    _config["mut-som-sig-mix"] = node;
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
    fprintf(stderr, "################################################################################\n");
    fprintf(stderr, " / ___| | ___  _ __ (_)  _ \\| |__  _   _ \n");
    fprintf(stderr, "| |   | |/ _ \\| '_ \\| | |_) | '_ \\| | | |\n");
    fprintf(stderr, "| |___| | (_) | | | | |  __/| | | | |_| |\n");
    fprintf(stderr, " \\____|_|\\___/|_| |_|_|_|   |_| |_|\\__, |\n");
    fprintf(stderr, "                                   |___/  (%s)\n", version::GIT_TAG_NAME);
    //fprintf(stderr, "%s %s\n", PROGRAM_NAME, version::GIT_TAG_NAME);
    fprintf(stderr, "================================================================================\n");
    fprintf(stderr, "Running with the following options:\n");
    fprintf(stderr, "================================================================================\n");
    fprintf(stderr, "  random seed:\t\t%ld\n", seed);
    if (fn_tree.length()>0) {
      fprintf(stderr, "  clone tree:\t\t%s\n", fn_tree.c_str());
    } else {
      fprintf(stderr, "  clones:\t\t%d\n", n_clones);
    }
    fprintf(stderr, "--------------------------------------------------------------------------------\n");
    fprintf(stderr, "Reference genome:\n");
    fprintf(stderr, "--------------------------------------------------------------------------------\n");
    if (fn_ref_fa.length()>0) {
      fprintf(stderr, "  input file:\t\t%s\n", fn_ref_fa.c_str());
    } else {
      fprintf(stderr, "  generate in-silico:\tyes\n");
      fprintf(stderr, "  ref seqs:\t\t%d\n", this->getValue<int>("ref-seq-num"));
      fprintf(stderr, "  seq len:\t\t%d (+/-%d)\n", this->getValue<int>("ref-seq-len-mean"), this->getValue<int>("ref-seq-len-mean"));
    }
    fprintf(stderr, "--------------------------------------------------------------------------------\n");
    fprintf(stderr, "Germline mutations:\n");
    fprintf(stderr, "--------------------------------------------------------------------------------\n");
    if (fn_mut_gl_vcf.length() > 0) {
      fprintf(stderr, "  germline VCF:\t%s\n", fn_mut_gl_vcf.c_str());
    } else {
      // TODO: print germline mutation params
      fprintf(stderr, "  number of mutations: %d\n", n_mut_gl);
      fprintf(stderr, "  evolutionary model:\t%s\n", mut_gl_model.c_str());
    }
    fprintf(stderr, "--------------------------------------------------------------------------------\n");
    fprintf(stderr, "Sequencing data\n");
    fprintf(stderr, "--------------------------------------------------------------------------------\n");
    if (fn_bam_input.length()>0) {
      fprintf(stderr, "  importing reads from BAM FILE provided by user\n");
      fprintf(stderr, "  input reads:\t%s\n", fn_bam_input.c_str());
    } else if (!seq_read_gen) {
      fprintf(stderr, "  simulating READ COUNTS\n");
      fprintf(stderr, "  seq depth:\t\t%d\n", this->getValue<int>("seq-coverage"));
      fprintf(stderr, "  depth dispersion:\t%.1f\n", this->getValue<double>("seq-rc-disp"));
      fprintf(stderr, "  seq error:\t\t%.2f\n", this->getValue<double>("seq-rc-error"));
      fprintf(stderr, "  min ALT read count:\t%d\n", this->getValue<int>("seq-rc-min"));
    } else {
      fprintf(stderr, "  simulating SEQUENCING READS\n");
      fprintf(stderr, "  reuse healthy reads:\t%s\n", do_reuse_reads ? "yes" : "no");
      fprintf(stderr, "  seq depth:\t\t%d\n", this->getValue<int>("seq-coverage"));
      fprintf(stderr, "  seq read length:\t%d\n", this->getValue<int>("seq-read-len"));
      fprintf(stderr, "  seq insert size:\t%d (+-%d)\n", this->getValue<int>("seq-frag-len-mean"), this->getValue<int>("seq-frag-len-sd"));
      fprintf(stderr, "  simulator:\t\t%s\n", this->getValue<string>("seq-art-path").c_str());
    }
    fprintf(stderr, "--------------------------------------------------------------------------------\n");
    fprintf(stderr, "Somatic mutations\n");
    fprintf(stderr, "--------------------------------------------------------------------------------\n");
    fprintf(stderr, "  number of mutations:\t%d\n", n_mut_somatic);
    fprintf(stderr, "  trunk mutations:\t%d\n", n_mut_trunk);
    if (_config["samples"]) {
      fprintf(stderr, "--------------------------------------------------------------------------------\n");
      fprintf(stderr, "Sampling scheme:\n");
      fprintf(stderr, "--------------------------------------------------------------------------------\n");
      map<string, vector<double>> sample_mtx = this->getMatrix<double>("samples");
      string s_sampling = stringio::printMatrix(sample_mtx);
      fprintf(stderr, "%s", s_sampling.c_str());
    }
    fprintf(stderr, "################################################################################\n");
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
