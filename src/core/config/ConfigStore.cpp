#include "ConfigStore.hpp"
#include "DataFrame.hpp"
#include "../stringio.hpp"

using namespace std;
using model::DataFrame;
namespace fs = boost::filesystem;

namespace config {

SampleConfig::SampleConfig(string label, vector<double> row) 
: m_label(label), m_vec_prevalence(row)
{}

// default constructor
ConfigStore::ConfigStore()
{
  // default values
  m_mut_gl_model = "JC";
  m_mut_gl_model_params_kappa = 1.0;
  m_mut_gl_model_params_nucfreq = vector<double>(4, 0.25);

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
  string dir_out = "output";
  string fn_config = "";
  string fn_bam_input = "";
  string fn_mut_gl_vcf = "";
  string fn_mut_som_vcf = "";
  string fn_mut_som_sig = "resources/signatures_probabilities.txt";
  string fn_ref_trinuc_sig = "";
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

  // program description
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
    ("ref-fasta,r", po::value<string>(&fn_ref_fa), "reference sequence")
    ("mut-gl-vcf,v", po::value<string>(&fn_mut_gl_vcf), "germline variants")
    ("mut-som-trunk,i", po::value<int>(&n_mut_trunk), "number of transforming mutations (separating healthy genome from first cancer genome)")
    ("tree,t", po::value<string>(&fn_tree), "file containing user defined clone tree (Newick format)")
    ("out-dir,o", po::value<string>(&dir_out), "output directory")
    ("verbosity,v", po::value<int>(&verb), "detail level of console output")
    ("threads,p", po::value<int>(&threads)->default_value(1), "number of parallel threads")
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
    fn_config = var_map["config"].as<string>();
    if (!fileExists(fn_config)) {
      fprintf(stderr, "\nArgumentError: File '%s' does not exist.\n", fn_config.c_str());
      return false;
    }
    // initialize global configuration from config file
    _config = YAML::LoadFile(fn_config);
  }

  // Manage paths / filenames
  //---------------------------------------------------------------------------

  // get executable's absolute directory
  fs::path path_exec( fs::initial_path<fs::path>() );
  path_exec = fs::system_complete( fs::path( av[0] ) ).parent_path();
  // get working directory
  fs::path path_work( fs::current_path() );
  // get config file's absolute directory
  fs::path path_conf( fs::initial_path<fs::path>() );
  if ( fn_config.length() > 0 ) { // find files relative to config directory
    path_conf = fs::path( fn_config );
    if ( path_conf.is_relative() ) {
      path_conf = fs::system_complete( path_conf );
    }
    path_conf = path_conf.parent_path();
  } else { // find files relative to current working directory
    path_conf = path_work;
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

  //---------------------------------------------------------------------------
  // reference-related params
  //---------------------------------------------------------------------------
  
  // path to reference FASTA
  if (var_map.count("ref-fasta") || !_config["ref-fasta"]) {
    _config["ref-fasta"] = fn_ref_fa;
  }
  fn_ref_fa = _config["ref-fasta"].as<string>();
  // path to reference trinucleotide profile
  if (var_map.count("ref-trinuc-profile") || !_config["ref-trinuc-profile"]) {
    _config["ref-trinuc-profile"] = fn_ref_trinuc_sig;
  }
  fn_ref_trinuc_sig = _config["ref-trinuc-profile"].as<string>();

  //---------------------------------------------------------------------------
  // mutation-related params
  //---------------------------------------------------------------------------

  // path to somatic mutation signature file
  if (!_config["mut-som-sig-file"]) {
    _config["mut-som-sig-file"] = fn_mut_som_sig;
  } 
  fn_mut_som_sig = _config["mut-som-sig-file"].as<string>();
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
  // sampling-related params
  //---------------------------------------------------------------------------

  // if sampling scheme was not provided, bail out
  if ( !_config["sampling"] || _config["sampling"].as<string>().size() == 0 ) {
    fprintf(stderr, "\nArgumentError: Parameter 'sampling' is required.\n");
    return false;
  } 
  // check if file name was provided
  else if ( _config["sampling"].Type() == YAML::NodeType::Scalar ) {
    fs::path path_prev( _config["sampling"].as<string>() );
    if ( path_prev.is_relative() ) { // relative path: find from config location
      path_prev = path_conf / path_prev;
    }
    // make sure file exists
    if ( !fileExists(path_prev.string()) ) {
      fprintf(stderr, "\nArgumentError: File '%s' does not exist.\n", path_prev.c_str());
      return false;
    }
    // read sampling scheme from CSV file
    vector<vector<string>> mtx_sampling;
    stringio::readCSV(mtx_sampling, path_prev.string());
    this->df_sampling = DataFrame<double>(mtx_sampling);
  } 
  // check if sampling matrix has been provided in config file
  else if (_config["sampling"].Type() == YAML::NodeType::Sequence) {
    // read sampling scheme from YAML node
    vector<vector<string>> mtx_sampling = getMatrix<string>("sampling");
    this->df_sampling = DataFrame<double>(mtx_sampling);
  }

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
  fs::path p_mut_som_sig( fn_mut_som_sig );
  if ( p_mut_som_sig.is_relative() ) { // relative path: find from executable location
    p_mut_som_sig = path_exec / p_mut_som_sig;
  }
  if ( fn_mut_som_sig.length()>0 && !fileExists( p_mut_som_sig.string() ) ) {
    fprintf(stderr, "\nArgumentError: Somatic mutation signature file '%s' does not exist.\n", p_mut_som_sig.c_str());
    return false;
  }
  _config["mut-som-sig-file"] = p_mut_som_sig.string();
  // was a clone tree provided by the user?
  fs::path path_tree( fn_tree );
  if ( fn_tree.length()>0 ) {
    if ( path_tree.is_relative() ) { // relative path: find from config location
      path_tree = path_conf / path_tree;
    }
    // check: does tree file exist?
    if ( !fileExists( path_tree.string() ) ) {
      fprintf(stderr, "\nArgumentError: Tree file '%s' does not exist.\n", path_tree.c_str());
      return false;
    }
  }
  _config["tree"] = path_tree.string();

  // check germline evolutionary model params
  if (!_config["mut-gl-model"] || _config["mut-gl-model"].as<string>().size() == 0) {
    //fprintf(stderr, "\n[INFO] Missing evolutionary model - assuming '%s'.\n", model.c_str());
    YAML::Node node = YAML::Load(m_mut_gl_model);
    _config["mut-gl-model"] = node;
  }
  m_mut_gl_model = _config["mut-gl-model"].as<string>();
  if (!_config["mut-gl-model-params"]) {
    //fprintf(stderr, "\n[INFO] Missing evolutionary model - assuming '%s'.\n", model.c_str());
    YAML::Node node;
    node["kappa"] = m_mut_gl_model_params_kappa;
    node["nucFreq"] = m_mut_gl_model_params_nucfreq;
    _config["mut-gl-model-params"] = node;
  }
  m_mut_gl_model_params_kappa   = _config["mut-gl-model-params"]["kappa"].as<double>();
  m_mut_gl_model_params_nucfreq = _config["mut-gl-model-params"]["nucFreq"].as<vector<double>>();
  if (m_mut_gl_model != "JC" && !_config["mut-gl-model-params"]) {
    fprintf(stderr, "\nArgumentError: Evolutionary models other than 'JC' require parameters (set param 'mut-gl-model-params' in config file.)\n");
    return false;
  }
  if (m_mut_gl_model == "F81" || m_mut_gl_model == "HKY") {
    if (!_config["mut-gl-model-params"]["nucFreq"]) {
      fprintf(stderr, "\nArgumentError: Model '%s' requires nucleotide frequencies (set param 'mut-gl-model-params':'nucFreq' in config file.)\n", m_mut_gl_model.c_str());
      return false;
    } else if (_config["mut-gl-model-params"]["nucFreq"].size() != 4) {
      fprintf(stderr, "\nArgumentError: Parameter 'mut-gl-model-params':'nucFreq' must contain exactly 4 values.\n");
      return false;
    }
  }
  if (m_mut_gl_model == "K80" || m_mut_gl_model == "HKY") {
    if (!_config["mut-gl-model-params"]["kappa"]) {
      fprintf(stderr, "\nArgumentError: Model '%s' requires transition-transversion ratio (set param 'mut-gl-model-params':'kappa' in config file.)\n", m_mut_gl_model.c_str());
      return false;
    }
  }
  if (m_mut_gl_model == "matrix") {
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
  if (!_config["mut-som-sig-mix"]) {
    fprintf(stderr, "\n[INFO] Missing somatic mutation signatures - assuming signature 1.\n");
    YAML::Node node = YAML::Load("{'Signature 1': 1.0}");
    _config["mut-som-sig-mix"] = node;
  }

  // does sampling matrix have the expected number of rows (clones + 1)?
  double row_sum;
  for (unsigned i=0; i<this->df_sampling.n_rows; i++) {
    // check: row sum <= 1?
    row_sum = 0.0;
    string row_label = this->df_sampling.rownames[i];
    vector<double> row = this->df_sampling.data[i];
    for (size_t j=1; j<row.size(); j++)
      row_sum += row[j];
    if (row_sum > 1) {
      fprintf(stderr, "\n[WARN]: Sampling matrix contains row sums >1 (sample '%s'). Normalizing...\n", row_label.c_str());
      for (size_t j=1; j<row.size(); j++) {
        row[j] /= row_sum;
      }
    }
    this->m_vec_samples.push_back(SampleConfig(row_label, row));
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
      fprintf(stderr, "  evolutionary model:\t%s\n", m_mut_gl_model.c_str());
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
    if (_config["sampling"]) {
      fprintf(stderr, "--------------------------------------------------------------------------------\n");
      fprintf(stderr, "Sampling scheme:\n");
      fprintf(stderr, "--------------------------------------------------------------------------------\n");
      fprintf(stderr, "%s", df_sampling.to_string().c_str());
    }
    fprintf(stderr, "################################################################################\n");
  }

  return true;
}

DataFrame<double> 
ConfigStore::getSamplingScheme ()
const 
{
  return this->df_sampling;
}

bool fileExists(string filename) {
  struct stat buffer;
  if (stat(filename.c_str(), &buffer)!=0) {
    return false;
  }
  return true;
}

} /* namespace config */
