/**
 * Simulation of clonal DNA sequences based on the process of somatic evolution.
 *
 * @author: Harald Detering (harald.detering@gmail.com)
 *
 */
#include "core/bamio.hpp"
#include "core/bamio/BulkSampleGenerator.hpp"
#include "core/clone.hpp"
#include "core/config/ConfigStore.hpp"
#include "core/model/DataFrame.hpp"
#include "core/random.hpp"
#include "core/seqio.hpp"
#include "core/treeio.hpp"
#include "core/vario.hpp"

#include "pcg-cpp/pcg_random.hpp"

#include <boost/filesystem.hpp> // absolute(),
#include <algorithm> // find()
#include <cassert>
#include <cstdio> // remove(), rename()
#include <ctime>
#include <exception>
#include <fstream>
#include <iostream>
#include <iterator> // begin(), end()
#include <map>
#include <math.h>
#include <memory>
#include <random> // std::random_device
#include <sstream>
#include <stdlib.h> // system()
#include <string>

using namespace std;
using namespace boost::filesystem;
using config::ConfigStore;
using model::DataFrame;
using stringio::format;
using seqio::SeqRecord;
using seqio::GenomeReference;
using seqio::GenomeInstance;
using vario::Genotype;
using vario::Mutation;
using vario::Variant;
using vario::VariantSet;
using vario::VariantStore;
using evolution::GermlineSubstitutionModel;
using evolution::SomaticSubstitutionModel;
using evolution::SomaticCnvModel;

// global constants
const string lbl_clone_normal = "N"; // TODO: make this a config param?

int main (int argc, char* argv[])
{
  // user params (defined in config file or command line)
  ConfigStore config;
  bool args_ok = config.parseArgs(argc, argv);
  if (!args_ok) { return EXIT_FAILURE; }

  // params specified by user (in config file or command line)
  path fn_bam_input(config.getValue<string>("bam-input"));
  int n_clones = config.getValue<int>("clones");
  int n_mut_somatic = config.getValue<int>("mut-som-num");
  int n_mut_transforming = config.getValue<int>("mut-som-trunk");
  int n_mut_germline = config.getValue<int>("mut-gl-num");
  double mut_gl_hom = config.getValue<double>("mut-gl-hom");
  unsigned ref_seq_num = config.getValue<int>("ref-seq-num");
  unsigned long ref_seq_len_mean = config.getValue<unsigned long>("ref-seq-len-mean");
  unsigned long ref_seq_len_sd = config.getValue<unsigned long>("ref-seq-len-sd");
  vector<double> ref_nuc_freqs = config.getValue<vector<double>>("ref-nuc-freqs");
  string fn_ref_fa = config.getValue<string>("ref-fasta");
  string fn_ref_trinuc_sig = config.getValue<string>("ref-trinuc-profile");
  string fn_mut_gl_vcf = config.getValue<string>("mut-gl-vcf");
  string str_model_gl = config.m_mut_gl_model;
  string fn_tree = config.getValue<string>("tree");
  string fn_mut_som_vcf = config.getValue<string>("mut-som-vcf");
  string fn_mut_som_sig = config.getValue<string>("mut-som-sig-file");
  //map<string, vector<double>> mtx_sample = config.getMatrix<double>("samples");
  DataFrame<double> df_sampling = config.getSamplingScheme();
  double seq_coverage = config.getValue<double>("seq-coverage");
  double seq_rc_error = config.getValue<double>("seq-rc-error");
  double seq_rc_disp = config.getValue<double>("seq-rc-disp");
  int    seq_rc_min = config.getValue<int>("seq-rc-min");
  bool seq_read_gen = config.getValue<bool>("seq-read-gen");
  bool seq_use_vaf = config.getValue<bool>("seq-use-vaf");
  unsigned seq_read_len = config.getValue<unsigned>("seq-read-len");
  unsigned seq_frag_len_mean = config.getValue<unsigned>("seq-frag-len-mean");
  unsigned seq_frag_len_sd = config.getValue<unsigned>("seq-frag-len-sd");
  unsigned seq_tile_pad = 0; //seq_frag_len_mean;
  string seq_art_path = config.getValue<string>("seq-art-path");
  bool seq_reuse_reads = config.getValue<bool>("seq-reuse-reads");
  bool seq_fq_out = config.getValue<bool>("seq-fq-out");
  bool seq_sam_out = config.getValue<bool>("seq-sam-out");
  int verbosity = config.getValue<int>("verbosity");
  int num_threads = config.threads;
  long seed = config.getValue<long>("seed");
  //string pfx_out = config.getValue<string>("out-pfx");
  string dir_out = config.getValue<string>("out-dir");
  double mut_som_cnv_ratio = config.getValue<double>("mut-som-cnv-ratio");

  //float ado_pct = 0.0; // TODO: make this a user parameter
  //int ado_frag_len = 10000; // TODO: make this a user parameter

  // inferred properties
  bool do_ref_sim = ( fn_ref_fa.length() == 0 );
  bool do_ref_trinuc_sig = ( fn_ref_trinuc_sig.length() > 0 );
  bool do_germline_vars = ( (fn_mut_gl_vcf.size() > 0) || (n_mut_germline > 0) );
  bool do_somatic_vars = ( fn_mut_som_vcf.size() == 0 );
  bool do_cnv_sim = ( mut_som_cnv_ratio > 0.0 );
  // has the user provided a BAM file?
  bool seq_user_bam = fn_bam_input.string().length() == 0;
  // when simulating CNVs, reads cannot be reused (varying coverage required)
  seq_reuse_reads = !do_cnv_sim && seq_reuse_reads;
  // optimization: ref genome reads are generated as baseline for samples
  bool do_ref_reads = !seq_user_bam && seq_reuse_reads;

  // internal vars
  map<shared_ptr<Clone>, string> clone2fn; // stores genome file names for each clone
  treeio::Tree<Clone> tree; // contains clone genealogy
  GermlineSubstitutionModel model_gl; // model of germline seq evolution
  SomaticSubstitutionModel model_sm; // model of somatic seq evolution
  vector<Mutation> vec_mut_gl; // germline mutations
  vector<Mutation> vec_mut_som = vector<Mutation>(n_mut_somatic); // somatic mutations
  vector<Variant> vec_var_gl; // germline variants
  VariantSet varset_gl; // germline variants
  VariantSet varset_sm; // somatic variants
  // TODO: check previous variables, maybe consolidate functionality?
  VariantStore var_store; // manages somatic variants
  bamio::BulkSampleGenerator bulk_generator; // simulates bulk samples

  // set number of parallel threads
  omp_set_num_threads(num_threads);

  // inititalize germline model of sequence evolution
  if (str_model_gl == "JC") {
    model_gl.init_JC();
  } else if (str_model_gl == "F81") {
    //vector<double> vec_freqs = config.getValue<vector<double>>("mut-gl-model-params:nucFreq");
    vector<double> vec_freqs = config.m_mut_gl_model_params_nucfreq;
    model_gl.init_F81(&vec_freqs[0]);
  } else if (str_model_gl == "K80") {
    //double kappa = config.getValue<double>("mut-gl-model-params:kappa");
    double kappa = config.m_mut_gl_model_params_kappa;
    model_gl.init_K80(kappa);
  } else if (str_model_gl == "HKY") {
    //string kappa = config.getValue<string>("mut-gl-model-params:kappa");
    //vector<double> vec_freqs = config.getValue<vector<double>>("mut-gl-model-params:nucFreq");
    double kappa = config.m_mut_gl_model_params_kappa;
    vector<double> vec_freqs = config.m_mut_gl_model_params_nucfreq;
    double pi_i[4];
    copy(vec_freqs.begin(), vec_freqs.end(), pi_i);
    model_gl.init_HKY(pi_i, kappa);
  }
  // print substitution rate matrix for germline mutations?
  if (verbosity >= 2) {
    map<string, vector<double>> Q;
    int i = 0;
    for (auto row = begin(model_gl.Q); row!=end(model_gl.Q); ++row, i++) {
      //string nuc_from = stringio::format("%s", seqio::idx2nuc(i));
      stringstream ss;
      string nuc_from;
      ss << seqio::idx2nuc(i);
      ss >> nuc_from;
      Q[nuc_from] = vector<double>(begin(*row), end(*row));
    }
    fprintf(stdout, "Germline nucleotide substitution rates:\n%s", stringio::printMatrix(Q).c_str());
  }

  // initialize somatic model of sequence evolution
  if (do_somatic_vars) {
    map<string, double> map_mut_som_sig_mix = config.getMap<double>("mut-som-sig-mix");
    model_sm = SomaticSubstitutionModel(fn_mut_som_sig, map_mut_som_sig_mix);
  }
  SomaticCnvModel model_cnv; // model of somatic structural evolution
  if (do_cnv_sim) {
    model_cnv.len_min = config.getValue<double>("mut-som-cnv-len-min");
    model_cnv.len_exp = config.getValue<double>("mut-som-cnv-len-exp");
    model_cnv.gain_prob = config.getValue<double>("mut-som-cnv-gain-prob");
    model_cnv.rate_wgd = config.getValue<double>("mut-som-cnv-rate-wgd");
    model_cnv.rate_chr = config.getValue<double>("mut-som-cnv-rate-chr");
    model_cnv.rate_arm = config.getValue<double>("mut-som-cnv-rate-arm");
    model_cnv.rate_tel = config.getValue<double>("mut-som-cnv-rate-tel");
    model_cnv.rate_foc = config.getValue<double>("mut-som-cnv-rate-foc");
  }

  // create output directories
  path path_out(dir_out);
  path_out = absolute(path_out);
  path path_fasta = path_out / "fasta";
  path path_bam = path_out / "bam";
  path path_bed = path_out / "bed";
  if ( ! exists(path_out) ) {
    if ( ! create_directories(path_out) ) {
      cerr << "[ERROR] (main) Could not create output directory: " << path_out << endl;
      cerr << "[ERROR] (main) Bailing out..." << endl;
      return EXIT_FAILURE;
    }
    create_directories(path_fasta);
    create_directories(path_bam);
    create_directories(path_bed);
  } else { // output dir exists
    // bail out with error
    //cerr << "[ERROR] (main) Output path already exists: " << path_out << "" << endl;
    //cerr << "[ERROR] (main) Please change/remove it and restart." << endl;
    //return EXIT_FAILURE;

    cerr << "[WARN] (main) Output path already exists: " << path_out << "" << endl;

    // create subdirectories if they don't exist
    if ( ! exists(path_fasta) ) {
      create_directories(path_fasta);
    }
    if ( ! exists(path_bam) ) {
      create_directories(path_bam);
    }
    if ( ! exists(path_bed) ) {
      create_directories(path_bed);
    }
  }

  // output file names
  const path fn_mm = path_out / "mm.csv"; // mutation map
  const path fn_newick = path_out / "clone.tree"; // clone tree
  const path fn_dot = path_out / "clone.tree.dot"; // DOT graph for clone tree
  const path fn_sampling_csv = path_out / "prevalences.csv"; // sampling scheme

  // initialize random functions
  //seed = time(NULL) + clock();
  fprintf(stderr, "random seed: %ld\n", seed);
  //RandomNumberGenerator<> rng(seed);
  RandomNumberGenerator rng(seed);
  function<double()> random_dbl = rng.getRandomFunctionReal(0.0, 1.0);
  function<double()> random_gamma = rng.getRandomFunctionGamma(2.0, 0.25);

  if (fn_tree.size() > 0) {
    cerr << "Reading tree from file '" << fn_tree << "'" << endl;
    tree = *(new treeio::Tree<Clone>(fn_tree));
    n_clones = tree.getVisibleNodes().size();
    cerr << "num_nodes:\t" << tree.m_numVisibleNodes << endl;
    cerr << "n_clones:\t" << n_clones << endl;

    // if number of mutations have not been supplied specifically,
    // branch lengths are interpreted as number of mutations
    if (n_mut_somatic == 0) {
      fprintf(stderr, "\nNumber of mutations not specified, using tree branch lengths (expected number of mutations).\n");
      // TODO: should absolute number of mutations be encodable in branch lengths?
      double dbl_branch_length = tree.getTotalBranchLength();
      n_mut_somatic = floor(dbl_branch_length);
      fprintf(stderr, "Simulating a total of %d mutations.\n", n_mut_somatic);
    }
  }
  else if (n_clones > 0) {
    tree = treeio::Tree<Clone>(n_clones);
    fprintf(stderr, "\nGenerating random topology for %d clones...\n", n_clones);
    tree.generateRandomTopologyLeafsOnly(random_dbl);
    tree.varyBranchLengths(random_gamma);
    fprintf(stderr, "done.\n");
  }

  fprintf(stderr, "\nNewick representation of underlying tree:\n");
  tree.printNewick(cerr);
  fprintf(stderr, "\n");

  // drop somatic mutations on clone tree
  tree.dropSomaticMutations(n_mut_somatic, n_mut_transforming, rng);

  fprintf(stderr, "\nWriting mutated tree to file (Newick format): %s\n", fn_newick.c_str());
  tree.printNewick(fn_newick.string());

  fprintf(stderr, "Writing mutated tree to file (DOT format): %s\n", fn_dot.c_str());
  tree.printDot(fn_dot.string());

  // TODO: should output contain only SNVs here?
  // write mutation matrix to file
  fprintf(stderr, "Writing mutation matrix for visible clones to file: %s\n", fn_mm.c_str());
  tree.writeMutationMatrix(fn_mm.string());


  // get reference genome
  GenomeReference ref_genome;
  if (fn_ref_fa.length() > 0) {
    fprintf(stderr, "\nReading reference from file '%s'...\n", fn_ref_fa.c_str());
    ref_genome = GenomeReference(fn_ref_fa.c_str());
  }
  else if (ref_nuc_freqs.size() > 0) {
    fprintf(stderr, "\nGenerating random reference genome sequence (%u seqs of length %lu +/- %lu)...", ref_seq_num, ref_seq_len_mean, ref_seq_len_sd);
    ref_genome.generate(ref_seq_num, ref_seq_len_mean, ref_seq_len_sd, ref_nuc_freqs, rng);

    // write generated genome to file
    fn_ref_fa = (path_out / path("ref.fa")).string();
    fprintf(stderr, "writing reference genome to file: %s\n", fn_ref_fa.c_str());
    seqio::writeFasta(ref_genome.records, fn_ref_fa);
    clone2fn[tree.m_root] = fn_ref_fa;
  }
  ref_genome.indexRecords();
  fprintf(stderr, "read (%u bp in %u sequences).\n", ref_genome.length, ref_genome.num_records);
  // TODO: make this a parameter?
  ref_genome.ploidy = 2;

  path fn_ref_bam;
  if (do_ref_reads) {
    // generate baseline sequencing reads from reference genome (haploid)
    fprintf(stderr, "---\nGenerating baseline sequencing reads from reference genome...\n");
    // initialize ART wrapper
    bamio::ArtWrapper art(seq_art_path);
    art.read_len = seq_read_len;
    art.frag_len_mean = seq_frag_len_mean;
    art.frag_len_sd = seq_frag_len_sd;
    art.fold_cvg = seq_coverage;
    art.fn_ref_fa = fn_ref_fa;
    art.out_sam = seq_sam_out;
    //art.do_keep_fq = seq_fq_out; // do baseline FASTQs have any application?
    // generate reads using ART
    path art_out_pfx = path_out / path("ref");
    int res_art = art.run(art_out_pfx.string());
    // Simulated read alignments can be found in the following SAM file
    fn_ref_bam = art_out_pfx;
    fn_ref_bam += ".sam";

    // simulate individual baseline reads for each sample if required
    // NOTE: if simulating CNVs, samples are created later from clone-specific read sets
    if (!seq_reuse_reads && !do_cnv_sim) {
      fprintf(stderr, "---\nGenerating baseline sequencing reads for %us samples...\n", df_sampling.n_rows);
      for (unsigned i=0; i<df_sampling.n_rows; i++) {
        string id_sample = df_sampling.rownames[i];
        art_out_pfx = path_out / path(format("%s.baseline", id_sample.c_str()));
        res_art = art.run(art_out_pfx.string());
      }
    }
  } else {
    fn_ref_bam = fn_bam_input;
  }

  // if a VCF file was provided, read germline variants from file, otherwise simulate variants
  if (fn_mut_gl_vcf.size() > 0) { // TODO: check consistency VCF <-> reference
    fprintf(stderr, "applying germline variants (from %s).\n", fn_mut_gl_vcf.c_str());
    VariantSet ref_variants;
    map<string, vector<Genotype >> ref_gt_matrix;
    vario::readVcf(fn_mut_gl_vcf, varset_gl, ref_gt_matrix);
  } else if (n_mut_germline > 0) {
    fprintf(stderr, "simulating %d germline variants (model: %s).\n", n_mut_germline, str_model_gl.c_str());
    //vec_var_gl = var_store.generateGermlineVariants(n_mut_germline, ref_genome, model_gl, rng);
    //varset_gl = VariantSet(vec_var_gl);
    var_store.generateGermlineVariants(n_mut_germline, ref_genome, model_gl, mut_gl_hom, rng);

    fn_mut_gl_vcf = format("%s/germline.vcf", path_out.c_str());
    fprintf(stderr, "writing generated variants to file: %s\n", fn_mut_gl_vcf.c_str());
    //vario::writeVcf(ref_genome.records, vec_var_gl, "germline", fn_mut_gl_vcf);
    var_store.writeGermlineSnvsToVcf(fn_mut_gl_vcf, ref_genome);
  }

  // assign somatic mutation types (single-nucleotide vs. copy-number)
  unsigned n_som_cnv = vario::assignSomaticMutationType(vec_mut_som, mut_som_cnv_ratio, rng);

  // generate somatic mutations
  vector<Variant> vec_var_somatic;
  if (do_somatic_vars) {
    // generate point mutations (relative position + chr copy)
    var_store.generateSomaticVariants(vec_mut_som, ref_genome, model_sm, model_cnv, rng);
    vec_var_somatic = var_store.getSomaticSnvVector();
    varset_sm = VariantSet(vec_var_somatic);

    // export variants to output files
    path fn_cnv = path_out / "somatic.cnv.tsv";
    var_store.writeCnvsToFile(fn_cnv.string());
//fprintf(stderr, "\nTotal set of mutations (id, rel_pos, copy):\n");
//for (int i=0; i<n_mut_somatic; i++)
//  fprintf(stderr, "%d\t%f\t%d\n", mutations[i].id, mutations[i].relPos, mutations[i].copy);
  }

  // get nodes from clone tree (in hierarchical order)
  vector<shared_ptr<Clone>> nodes = tree.getNodesPreOrder();
  // keep an individual GenomeInstance for each clone tree node
  map<int, GenomeInstance> map_id_genome;
  // generate "healthy" GenomeInstance from reference genome
  GenomeInstance healthy_genome = GenomeInstance(ref_genome);
  // root node corresponds to healthy genome (diploid)
  map_id_genome[nodes[0]->index] = healthy_genome;
  // apply germline mutations to healthy genome
  fprintf(stderr, "applying germline variants to clone tree root...\n");
  var_store.applyGermlineVariants(healthy_genome, rng);

  // generate GenomeInstances for clones
  if (nodes.size() > 1) {
    fprintf(stderr, "applying somatic variants by traversing clone tree...\n");

    for (size_t i=1; i<nodes.size(); ++i) {
      cerr << "\t" << *(nodes[i]) << endl;

      // copy parent's genome
      vector<seqio::seg_mod_t> vec_seg_mod;
      //GenomeInstance gi_node(map_id_genome[nodes[i]->parent->index], vec_seg_mod);
      GenomeInstance gi_parent = map_id_genome[nodes[i]->parent->index];
      GenomeInstance gi_node(gi_parent, vec_seg_mod);
//cerr << "parent:\t" << *(gi_parent.vec_chr[0]->lst_segments.begin()) << endl;
//cerr << "child:\t" << *(gi_node.vec_chr[0]->lst_segments.begin()) << endl;
      var_store.transferMutations(vec_seg_mod);

      // get mutations for clone tree branch
      vector<int> node_mut = nodes[i]->m_vec_mutations;
      fprintf(stderr, "\t\t[%lu mutations...]\n", node_mut.size());
      // apply somatic mutations (SNVs + CNVs, in order)
      for (int mut : node_mut) {
//cerr << "id_mut: " << mut << endl;
        var_store.applyMutation(vec_mut_som[mut], gi_node, rng);
      }
      // store GenomeInstance for further use
      map_id_genome[nodes[i]->index] = gi_node;
    }
  }

cerr << "parent:\t" << *(map_id_genome[0].vec_chr[0]->lst_segments.begin()) << endl;
cerr << "child:\t" << *(map_id_genome[1].vec_chr[0]->lst_segments.begin()) << endl;

  // prepare clone genomes for export
  map<string, GenomeInstance> map_clone_genome;
  vector<string> vec_clone_lbl;
  for (size_t i=1; i<nodes.size(); ++i) {
    if (!nodes[i]->is_visible) continue;
    int idx = nodes[i]->index;
    string lbl = nodes[i]->label;
    map_clone_genome[lbl] = map_id_genome[idx];
    vec_clone_lbl.push_back(lbl);
  }
  // sanity check: make sure labels in sampling scheme match clone labels
  for (auto lbl : df_sampling.colnames) {
    auto it = vec_clone_lbl.begin();
    while (it != vec_clone_lbl.end() && *it != lbl) it++;
    if (it == vec_clone_lbl.end()) {
      fprintf(stderr, "[ERROR] (main) Label '%s' not found in clone tree.\n", lbl.c_str());
      return EXIT_FAILURE;
    }
  }
  // add "healthy" clone
  shared_ptr<Clone> c_normal = nodes[0];
  map_clone_genome[lbl_clone_normal] = map_id_genome[c_normal->index];
  vec_clone_lbl.push_back(lbl_clone_normal);


// DEBUG info: export segment copies to file
// TODO: move out to Logger class
string fn_dbg_segs = (path_out / "segments.tsv").string();
std::ofstream ofs_dbg_segs(fn_dbg_segs, std::ofstream::out);
for (auto const & lbl_gi : map_clone_genome) {
  string lbl_clone = lbl_gi.first;
  for (auto const & id_chr : lbl_gi.second.map_id_chr) {
    string lbl_chr = id_chr.first;
    for (auto const & chr : id_chr.second) {
      for (auto const & seg : chr->lst_segments) {
        ofs_dbg_segs << lbl_clone << "\t";
        ofs_dbg_segs << lbl_chr << "\t";
        ofs_dbg_segs << seg.id << "\t";
        ofs_dbg_segs << seg.ref_start << "\t";
        ofs_dbg_segs << seg.ref_end << endl;
      }
    }    
  }
}
ofs_dbg_segs.close();

// DEBUG info: export variants associated to each segment copy
// TODO: move out to Logger class
string fn_dbg_vars = (path_out / "segment_vars.tsv").string();
std::ofstream ofs_dbg_vars(fn_dbg_vars, std::ofstream::out);
for (auto const & cg : map_clone_genome) {
  string lbl_clone = cg.first;
  for (auto const & ic : cg.second.map_id_chr) {
    string lbl_chr = ic.first;
    for (auto const & chr : ic.second) {
      for (auto const & seg : chr->lst_segments) {
        for (auto const & v : var_store.map_seg_vars[seg.id]) {
          ofs_dbg_vars << seg.id << "\t" << var_store.map_id_snv[v].id << endl;
        }
      }
    }
  }
}
ofs_dbg_vars.close();

  // for decoupling, transform sampling data frame to map
  map<string, map<string, double>> mtx_sample_clone;
  for (unsigned i=0; i<df_sampling.n_rows; i++) {
    string id_sample = df_sampling.rownames[i];
    vector<double> w = df_sampling.data[i];

    // sanity check: weights and clone labels must match
    assert ( w.size()+1 == vec_clone_lbl.size() );

    map<string, double> map_clone_weight;
    double sum_w = 0.0;
    for (size_t i=0; i<w.size(); i++) {
      map_clone_weight[df_sampling.colnames[i]] = w[i];
      sum_w += w[i];
    }

    // normalize clone weights?
    if (sum_w > 1.0) {
      for (size_t i=0; i<w.size(); i++) {
        map_clone_weight[df_sampling.colnames[i]] /= sum_w;
      }
      // set normal cell combination
      map_clone_weight[lbl_clone_normal] = 0.0;
    } 
    else { // no normalization, difference to 1 is normal contamination
      // set normal cell combination
      map_clone_weight[lbl_clone_normal] = 1.0 - sum_w;
    }

    //
    mtx_sample_clone[id_sample] = map_clone_weight;
  }

  // sanity check: number of genome matches clone fractions?
  assert ( map_clone_genome.size() == mtx_sample_clone.begin().second.size() );

  // initialize bulk samples
  bulk_generator.initSamples(mtx_sample_clone);

  // initialize BulkSampleGenerator with reference genome
  bulk_generator.initRefSeqs(ref_genome);

  // initialize BulkSampleGenerator with clone genomes
  // (copy number segmentation)
  bulk_generator.initCloneGenomes(
    map_clone_genome,
    path_bed
  );

  // export genomic sequences (only if reads are to be generated)
  if ( seq_read_gen ) {
    fprintf(stdout, "Writing tiled ref seqs...\n");
    bulk_generator.writeCloneGenomes (
      map_clone_genome,
      ref_genome, 
      seq_tile_pad,
      seq_frag_len_mean, 
      path_fasta, 
      path_bed);
  }

  // determine allele-specific copy number state for samples
  // bulk_generator.calculateBulkCopyNumber(mtx_sample_clone, map_clone_genome);
  bulk_generator.calculateBulkCopyNumber(map_clone_genome);
  // write absolute copy number states to BED file for each sample
  bulk_generator.writeBulkCopyNumber(path_bed);

  // generate reads for genomic regions
  // - reads overlapping with padded regions are discarded
  fprintf(stdout, "Generating sequencing reads...\n");
  bulk_generator.generateBulkSamples (
    mtx_sample_clone,
    var_store,
    path_fasta,
    path_bam,
    path_bed,
    seq_coverage,
    seq_rc_error,
    seq_rc_disp,
    seq_rc_min,
    seq_read_gen,
    seq_use_vaf,
    seq_read_len, 
    seq_frag_len_mean, 
    seq_frag_len_sd,  
    seq_art_path,
    rng
  );

  // setup mutation matrix
  int num_nodes = tree.m_numNodes;
  vector<vector<bool> > mat_mut(num_nodes, vector<bool>(n_mut_somatic, false));
  tree.m_root->populateMutationMatrixRec(mat_mut);

  // collect clone labels
  vector<int> vec_vis_nodes_idx = tree.getVisibleNodesIdx();
  vector<string> vec_labels;
  for (auto i : vec_vis_nodes_idx) {
    vec_labels.push_back(tree.m_vecNodes[i]->label);
  }
  // sanity check: do visible clone labels match prevalence matrix?
  set<string> set_labels_tree;
  set<string> set_labels_prev;

  // write somatic variants to VCF file
  string fn_vcf = (path_out / "somatic.vcf").string();
  fprintf(stderr, "\nWriting (sub)clonal mutations to file '%s'.\n", fn_vcf.c_str());
  std::ofstream f_vcf;
  f_vcf.open(fn_vcf);
  vario::writeVcf(ref_genome.records, vec_var_somatic, vec_vis_nodes_idx, vec_labels, mat_mut, f_vcf);
  f_vcf.close();

  // DEPRECATED: Now handled by BulkSampleGenerator :)
  // // get clones and mutation map from intermediate files
  // vector<shared_ptr<Clone>> clones = tree.getVisibleNodes();
  // int num_mm_rows = 0;
  // map<string, vector<bool>> mm;
  // num_mm_rows = vario::readMutMapFromCSV(mm, fn_mm.string());
  // assert( num_mm_rows == clones.size() );

  // // spike-in mutations in baseline reads
  // VariantSet varset_spikein;

  // // generate normal sample
  // path fn_normal_bam = "";
  // fprintf(stderr, "---\nGenerating sequencing reads for normal (healthy) sample...\n");
  // if (do_germline_vars) {
  //   varset_spikein = varset_gl + varset_sm;
  //   fn_normal_bam = path_out / "normal.sam";
  //   path fn_fqout = path_out / "normal.fastq";
  //   path fn_baseline = fn_ref_bam;
  //   string id_sample("normal");
  //   vector<double> vec_freq(n_clones-1, 0.0);
  //   fprintf(stderr, "\nSpike in variants from file: %s\n", fn_mut_gl_vcf.c_str());
  //   fprintf(stderr, "  into baseline reads: %s\n", fn_ref_bam.c_str());
  //   bamio::mutateReads(fn_fqout, fn_normal_bam, fn_baseline, varset_spikein,
  //     clones, mm, vec_freq, id_sample, ref_genome.ploidy, rng, seq_fq_out);
  //   fprintf(stderr, "\nReads containing germline mutations are in file: %s\n", fn_normal_bam.c_str());
  // } else { // there are no germline variants
  //   fn_normal_bam = fn_ref_bam;
  // }
  // fprintf(stderr, "---\nGenerating sequencing reads for %ld tumor samples...\n", mtx_sample.size());
  // for (auto sample : mtx_sample) {
  //   path fn_baseline;
  //   if (seq_reuse_reads) {
  //     // use normal sample reads as baseline for tumor samples
  //     fn_baseline = fn_normal_bam;
  //     // germline variants are already contained in baseline reads
  //     varset_spikein = varset_sm;
  //   }
  //   else  {
  //     // choose sample-specific baseline files
  //     fn_baseline = path_out / str(format("%s.baseline.sam") % sample.first);
  //     // include germline variants in spike-in variants
  //     varset_spikein = varset_gl + varset_sm;
  //   }
  //   path fn_fqout = path_out / str(format("%s.bulk.fq") % sample.first);
  //   path fn_samout = path_out / str(format("%s.bulk.sam") % sample.first);

  //   //VariantSet varset_spikein(vec_var_spikein);
  //   vector<shared_ptr<Clone>> vec_vis_clones = tree.getVisibleNodes();
  //   bamio::mutateReads(fn_fqout, fn_samout, fn_baseline, varset_spikein,
  //     clones, mm, sample.second, sample.first, ref_genome.ploidy, rng);
  // }

  return EXIT_SUCCESS;
}
