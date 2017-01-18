/**
 * Simulation of clonal DNA sequences based on the process of somatic evolution.
 *
 * @author: Harald Detering (harald.detering@gmail.com)
 *
 */
#include "core/bamio.hpp"
#include "core/clone.hpp"
#include "core/config/ConfigStore.hpp"
#include "core/random.hpp"
#include "core/seqio.hpp"
#include "core/treeio.hpp"
#include "core/vario.hpp"
#include "core/basicclonetree.hpp"
#include "core/coalescentclonetree.hpp"

#include <boost/format.hpp>
#include <cassert>
#include <ctime>
#include <exception>
#include <fstream>
#include <map>
#include <math.h>
#include <memory>
#include <sstream>
#include <stdlib.h> // system()
#include <string>

using namespace std;
using boost::format;
using boost::str;
using config::ConfigStore;
using seqio::SeqRecord;
using seqio::Genome;
using vario::Genotype;
using vario::Variant;
using vario::VariantSet;
using evolution::GermlineSubstitutionModel;
using evolution::SomaticSubstitutionModel;


int main (int argc, char* argv[])
{
  // user params (defined in config file or command line)
  ConfigStore config;
  bool args_ok = config.parseArgs(argc, argv);
  if (!args_ok) { return EXIT_FAILURE; }

  // params specified by user (in config file or command line)
  string fn_bam_input = config.getValue<string>("bam-input");
  int num_clones = config.getValue<int>("clones");
  int n_mut_clonal = config.getValue<int>("mutations");
  int n_mut_transforming = config.getValue<int>("init-muts");
  int n_mut_germline = config.getValue<int>("mut-germline");
  unsigned ref_seq_num = config.getValue<int>("ref-seq-num");
  unsigned long ref_seq_len_mean = config.getValue<unsigned long>("ref-seq-len-mean");
  unsigned long ref_seq_len_sd = config.getValue<unsigned long>("ref-seq-len-sd");
  vector<double> ref_nuc_freqs = config.getValue<vector<double>>("ref-nuc-freqs");
  string fn_ref_fa = config.getValue<string>("reference");
  string fn_mut_ref_vcf = config.getValue<string>("reference-vcf");
  string fn_mut_som_vcf = config.getValue<string>("mut-som-vcf");
  string str_model_gl = config.getValue<string>("model");
  string fn_tree = config.getValue<string>("tree");
  map<string, vector<double>> mtx_sample = config.getMatrix<double>("samples");
  double seq_coverage = config.getValue<double>("seq-coverage");
  unsigned seq_read_len = config.getValue<unsigned>("seq-read-len");
  unsigned seq_frag_len_mean = config.getValue<unsigned>("seq-frag-len-mean");
  unsigned seq_frag_len_sd = config.getValue<unsigned>("seq-frag-len-sd");
  string seq_art_path = config.getValue<string>("seq-art-path");
  bool seq_reuse_reads = config.getValue<bool>("seq-reuse-reads");
  int verbosity = config.getValue<int>("verbosity");
  long seed = config.getValue<long>("seed");
  string pfx_out = config.getValue<string>("out");
  bool do_cnv_sim = false; // indicates if CNVs will be considered

  float ado_pct = 0.0; // TODO: make this a user parameter
  int ado_frag_len = 10000; // TODO: make this a user parameter

  // inferred properties
  bool do_ref_sim = fn_bam_input.length() == 0;
  bool do_ref_reads = fn_bam_input.length() == 0;
  bool do_germline_vars = (fn_mut_ref_vcf.size() > 0) || (n_mut_germline > 0);

  // output file names
  string fn_mm = str(format("%s.mm.csv") % pfx_out); // mutation map

  // internal vars
  map<shared_ptr<Clone>, string> clone2fn; // stores genome file names for each clone
  treeio::Tree<Clone> tree; // contains clone genealogy
  GermlineSubstitutionModel model_gl; // model of germline seq evolution
  SomaticSubstitutionModel model_sm; // model of somatic seq evolution

  vector<Mutation> vec_mut_gl; // germline mutations
  vector<Variant> vec_var_gl; // germline variants

  // initialize random functions
  //seed = time(NULL) + clock();
  fprintf(stderr, "random seed: %ld\n", seed);
  RandomNumberGenerator<> rng(seed);
  function<double()> random_dbl = rng.getRandomFunctionDouble(0.0, 1.0);
  function<double()> random_gamma = rng.getRandomGamma(2.0, 0.25);

  if (num_clones > 0) {
    if (fn_tree.size()>0) {
      cerr << "Reading tree from file '" << fn_tree << "'" << endl;
      //treeio::parse::node root_node;
      //treeio::parse::readNewick(fn_tree, root_node);
      //tree = *(new treeio::Tree<Clone>(root_node));
      tree = *(new treeio::Tree<Clone>(fn_tree));
      num_clones = tree.getVisibleNodes().size();
      cerr << "num_nodes:\t" << tree.m_numVisibleNodes << endl;
      cerr << "num_clones:\t" << num_clones << endl;

      // if number of mutations have not been supplied specifically,
      // branch lengths are interpreted as number of mutations
      if (n_mut_clonal == 0) {
        fprintf(stderr, "\nNumber of mutations not specified, using tree branch lengths (expected number of mutations).\n");
        // TODO: should absolute number of mutations be encodable in branch lengths?
        double dbl_branch_length = tree.getTotalBranchLength();
        n_mut_clonal = floor(dbl_branch_length);
        fprintf(stderr, "Simulating a total of %d mutations.\n", n_mut_clonal);
      }
    }
    else {
      tree = treeio::Tree<Clone>(num_clones);
      fprintf(stderr, "\nGenerating random topology for %d clones...\n", num_clones);
      tree.generateRandomTopologyLeafsOnly(random_dbl);
      tree.varyBranchLengths(random_gamma);
      fprintf(stderr, "done.\n");
    }
    fprintf(stderr, "\nNewick representation of underlying tree:\n");
    tree.printNewick(cerr);
    fprintf(stderr, "\n");

    tree.evolve(n_mut_clonal, n_mut_transforming, rng);
    string fn_newick = str(format("%s.tree") % pfx_out);
    fprintf(stderr, "\nWriting mutated tree to file (Newick format): %s\n", fn_newick.c_str());
    tree.printNewick(fn_newick);

    string fn_dot = str(format("%s.tree.dot") % pfx_out);
    fprintf(stderr, "Writing mutated tree to file (DOT format): %s\n", fn_dot.c_str());
    tree.printDot(fn_dot);

    // write mutation matrix to file
    fprintf(stderr, "Writing mutation matrix for visible clones to file: %s\n", fn_mm.c_str());
    tree.writeMutationMatrix(fn_mm);
  }

  // get reference genome
  Genome ref_genome;
  if (fn_ref_fa.length() > 0) {
    fprintf(stderr, "\nReading reference from file '%s'...\n", fn_ref_fa.c_str());
    ref_genome = Genome(fn_ref_fa.c_str());
  }
  else if (ref_nuc_freqs.size() > 0) {
    fprintf(stderr, "\nGenerating random reference genome sequence (%u seqs of length %lu +/- %lu)...", ref_seq_num, ref_seq_len_mean, ref_seq_len_sd);
    ref_genome.generate(ref_seq_num, ref_seq_len_mean, ref_seq_len_sd, ref_nuc_freqs, rng);

    // write generated genome to file
    fn_ref_fa = str(format("%s.ref.fa") % pfx_out.c_str());
    fprintf(stderr, "writing reference genome to file: %s\n", fn_ref_fa.c_str());
    seqio::writeFasta(ref_genome.records, fn_ref_fa);
    clone2fn[tree.m_root] = fn_ref_fa;
  }
  ref_genome.indexRecords();
  fprintf(stderr, "read (%u bp in %u sequences).\n", ref_genome.length, ref_genome.num_records);
  // duplicate genome (all loci homozygous reference)
  fprintf(stderr, "duplicating reference sequence (haploid -> diploid).\n");
  ref_genome.duplicate();
  //ref_genome.indexRecords(); // working with haploid reference, no re-indexing necessary

  string fn_ref_bam = "";
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
    // generate reads using ART
    string art_out_pfx = pfx_out + ".ref";
    int res_art = art.run(art_out_pfx);
    // Simulated read alignments can be found in the following SAM file
    fn_ref_bam = art_out_pfx + ".sam";

    // simulate individual baseline reads for each sample if required
    if (!seq_reuse_reads) {
      fprintf(stderr, "---\nGenerating baseline sequencing reads for %ld samples...\n", mtx_sample.size());
      for (auto sample : mtx_sample) {
        art_out_pfx = str(format("%s.%s.baseline") % pfx_out % sample.first);
        res_art = art.run(art_out_pfx);
      }
    }
  } else {
    fn_ref_bam = fn_bam_input;
  }

  /* CNV simulation will be performed further down the line
  if (do_cnv_sim) { // need to generate baseline reads respecting copy number states
    for (auto sample : mtx_sample) {
      // TODO: will need to be updated to support CNVs (non-uniform coverage)
      art_cmd += str(format(" -o %s.%s_reads") % pfx_out % sample.first);
      fprintf(stderr, "---\nRunning ART with the following paramters:\n%s\n", art_cmd.c_str());
      int res_art = system(art_cmd.c_str());
      if (res_art != 0) {
        fprintf(stderr, "[ERROR] ART call had non-zero return value:\n        %s\n", art_cmd.c_str());
        return EXIT_FAILURE;
      }
    }
  } else { // no CNVs: joint baseline read set for all samples
    art_cmd += str(format(" -o %s_reads") % pfx_out);
    fprintf(stderr, "---\nRunning ART with the following paramters:\n%s\n", art_cmd.c_str());
    int res_art = system(art_cmd.c_str());
    if (res_art != 0) {
      fprintf(stderr, "[ERROR] ART call had non-zero return value:\n        %s\n", art_cmd.c_str());
      return EXIT_FAILURE;
    }
  }
  */

  // inititalize model of sequence evolution
  if (str_model_gl == "JC") {
    model_gl.init_JC();
  } else if (str_model_gl == "F81") {
    vector<double> vec_freqs = config.getValue<vector<double>>("model-params:nucFreq");
    model_gl.init_F81(&vec_freqs[0]);
  } else if (str_model_gl == "K80") {
    double kappa = config.getValue<double>("model-params:kappa");
    model_gl.init_K80(kappa);
  } else if (str_model_gl == "HKY") {
    vector<double> vec_freqs = config.getValue<vector<double>>("model-params:nucFreq");
    double kappa = config.getValue<double>("model-params:kappa");
    model_gl.init_HKY(&vec_freqs[0], kappa);
  }

  // if a VCF file was provided, read germline variants from file, otherwise simulate variants
  if (fn_mut_ref_vcf.size() > 0) { // TODO: check consistency VCF <-> reference
    fprintf(stderr, "applying germline variants (from %s).\n", fn_mut_ref_vcf.c_str());
  } else if (n_mut_germline > 0) {
    fprintf(stderr, "simulating %d germline variants (model: %s).\n", n_mut_germline, str_model_gl.c_str());
    //vec_mut_gl = vario::generateMutations(n_mut_germline, random_dbl);
    vec_var_gl = vario::generateGermlineVariants(n_mut_germline, ref_genome, model_gl, rng);
    //vario::applyVariants(ref_genome, vec_var_gl);

    fn_mut_ref_vcf = str(format("%s.germline.vcf") % pfx_out.c_str());
    fprintf(stderr, "writing generated variants to file: %s\n", fn_mut_ref_vcf.c_str());
    vario::writeVcf(ref_genome.records, vec_var_gl, "germline", fn_mut_ref_vcf);
  }

  // apply germline variants
  string fn_normal_bam = "";
  if (do_germline_vars) {
    fn_normal_bam = str(format("%s.normal.sam") % pfx_out);
    fprintf(stderr, "\nSpike in variants from file: %s\n", fn_mut_ref_vcf.c_str());
    fprintf(stderr, "  into baseline BAM: %s\n", fn_ref_bam.c_str());
    VariantSet ref_variants;
    map<string, vector<Genotype >> ref_gt_matrix;
    vario::readVcf(fn_mut_ref_vcf, ref_variants, ref_gt_matrix);
    bamio::mutateReads(fn_normal_bam, fn_ref_bam, ref_variants, ref_genome.ploidy, rng);
    fprintf(stderr, "\nReads containing germline mutations are in file: %s\n", fn_normal_bam.c_str());
    //vario::applyVariants(ref_genome, ref_variants.vec_variants, ref_gt_matrix["0"]);
  } else { // there are no germline variants
    fn_normal_bam = fn_ref_bam;
  }

  vector<Mutation> mutations;
  if (n_mut_clonal > 0) {
    // generate point mutations (relative position + chr copy)
    mutations = vario::generateMutations(n_mut_clonal, random_dbl);
    fprintf(stderr, "\nTotal set of mutations (id, rel_pos, copy):\n");
    for (int i=0; i<n_mut_clonal; i++)
      fprintf(stderr, "%d\t%f\t%d\n", mutations[i].id, mutations[i].relPos, mutations[i].copy);
  }
  // apply mutations, creating variants in the process
  vector<Variant> variants = vector<Variant>();
  vario::applyMutations(mutations, ref_genome, model_gl, random_dbl, variants);

  // setup mutation matrix
  int num_nodes = tree.m_numNodes;
  vector<vector<bool> > mat_mut(num_nodes, vector<bool>(n_mut_clonal, false));
  tree.m_root->populateMutationMatrixRec(mat_mut);

  // collect clone labels
  vector<int> vec_vis_nodes_idx = tree.getVisibleNodesIdx();
  vector<string> vec_labels;
  for (auto i : vec_vis_nodes_idx) {
    vec_labels.push_back(tree.m_vecNodes[i]->label);
  }

  // write clonal variants to VCF file
  string fn_vcf = str(format("%s.somatic.vcf") % pfx_out);
  fprintf(stderr, "\nWriting (sub)clonal mutations to file '%s'.\n", fn_vcf.c_str());
  ofstream f_vcf;
  f_vcf.open(fn_vcf);
  vario::writeVcf(ref_genome.records, variants, vec_vis_nodes_idx, vec_labels, mat_mut, f_vcf);
  f_vcf.close();

  // get clones and mutation map from intermediate files
  vector<shared_ptr<Clone>> clones = tree.getVisibleNodes();
  int num_mm_rows = 0;
  map<string, vector<bool>> mm;
  num_mm_rows = vario::readMutMapFromCSV(mm, fn_mm);
  assert( num_mm_rows == clones.size() );

  // introduce somatic mutations in baseline reads
  fprintf(stderr, "---\nNow generating sequencing reads for %ld samples...\n", mtx_sample.size());
  for (auto sample : mtx_sample) {
    string fn_baseline;
    if (seq_reuse_reads) // use normal sample reads as baseline for tumor samples
      fn_baseline = str(format("build/%s") % fn_normal_bam.c_str());
    else // choose sample-specific baseline files
      fn_baseline = str(format("build/%s.%s.baseline.sam") % pfx_out.c_str() % sample.first);
    string fn_fqout = str(format("%s.%s.bulk.fq") % pfx_out % sample.first);;
    string fn_samout = str(format("%s.%s.bulk.sam") % pfx_out % sample.first);;

    VariantSet varset_somatic(variants);
    vector<shared_ptr<Clone>> vec_vis_clones = tree.getVisibleNodes();
    bamio::mutateReads(fn_fqout, fn_samout, fn_baseline, varset_somatic,
      clones, mm, sample.second, sample.first, ref_genome.ploidy, rng);
  }

  // perform ADO -> on hold (better: generate coverage distribution for read sim)
  /*if (ado_pct > 0)
    for (auto kv : clone2fn)
      seqio::simulateADO(kv.second, ref_genome.masked_length, ado_pct, ado_frag_len, random_dbl);
  */

  return EXIT_SUCCESS;
}
