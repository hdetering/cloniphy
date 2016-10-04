/**
 * Simulation of clonal DNA sequences based on the process of somatic evolution.
 *
 * @author: Harald Detering (harald.detering@gmail.com)
 *
 */
#include "core/clone.hpp"
#include "core/config/ConfigStore.hpp"
#include "core/random.hpp"
#include "core/seqio.hpp"
#include "core/treeio.hpp"
#include "core/vario.hpp"
#include "core/basicclonetree.hpp"
#include "core/coalescentclonetree.hpp"

#include <boost/format.hpp>
#include <ctime>
#include <exception>
#include <fstream>
#include <map>
#include <math.h>
#include <memory>
#include <sstream>
#include <string>

using namespace std;
using config::ConfigStore;
using seqio::SeqRecord;
using seqio::Genome;
using vario::Genotype;
using vario::Variant;
using vario::VariantSet;
using evolution::SubstitutionModel;


int main (int argc, char* argv[])
{
  // user params (defined in config file or command line)
  ConfigStore config;
  bool args_ok = config.parseArgs(argc, argv);
  if (!args_ok) { return EXIT_FAILURE; }

  // params specified by user (in config file or command line)
  int num_clones = config.getValue<int>("clones");
  int num_mutations = config.getValue<int>("mutations");
  int num_transmuts = config.getValue<int>("init-muts");
  unsigned seq_num = config.getValue<int>("seq-num");
  unsigned long seq_len_mean = config.getValue<unsigned long>("seq-len-mean");
  unsigned long seq_len_sd = config.getValue<unsigned long>("seq-len-sd");
  vector<double> ref_nuc_freqs = config.getValue<vector<double>>("ref-nuc-freqs");
  string fn_ref = config.getValue<string>("reference");
  string ref_vcf = config.getValue<string>("reference-vcf");
  string str_model = config.getValue<string>("model");
  string fn_tree = config.getValue<string>("tree");
  map<string, vector<double>> sample_mtx = config.getMatrix<double>("samples");
  int verbosity = config.getValue<int>("verbosity");
  long seed = config.getValue<long>("seed");

  float ado_pct = 0.1; // TODO: make this a user parameter
  int ado_frag_len = 10000; // TODO: make this a user parameter

  // internal vars
  map<shared_ptr<Clone>, string> clone2fn; // stores genome file names for each clone
  treeio::Tree<Clone> tree; // contains clone genealogy

  // initialize random functions
  //seed = time(NULL) + clock();
  fprintf(stderr, "random seed: %ld\n", seed);
  RandomNumberGenerator<> rng(seed);
  function<double()> random_dbl = rng.getRandomFunctionDouble(0.0, 1.0);
  function<double()> random_gamma = rng.getRandomGamma(2.0, 0.25);

  if (num_clones > 0) {
    if (fn_tree.size()>0) {
      cerr << "Reading tree from file '" << fn_tree << "'" << endl;
      treeio::parse::node root_node;
      treeio::parse::readNewick(fn_tree, root_node);
      tree = *(new treeio::Tree<Clone>(root_node));
      num_clones = tree.getVisibleNodes().size();
      cerr << "num_nodes:\t" << tree.m_numVisibleNodes << endl;
      cerr << "num_clones:\t" << num_clones << endl;

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
      tree = treeio::Tree<Clone>(num_clones);
      fprintf(stderr, "\nGenerating random topology for %d clones...\n", num_clones);
      tree.generateRandomTopologyLeafsOnly(random_dbl);
      tree.varyBranchLengths(random_gamma);
      fprintf(stderr, "done.\n");
    }
    fprintf(stderr, "\nNewick representation of underlying tree:\n");
    tree.printNewick(cerr);
    fprintf(stderr, "\n");

    tree.evolve(num_mutations, num_transmuts, rng);
    fprintf(stderr, "\nNewick representation of mutated tree:\n");
    tree.printNewick(cerr);
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
  }

  // get reference genome
  Genome ref_genome;
  if (fn_ref.length() > 0) {
    fprintf(stderr, "\nGenerating random reference genome sequence (%u seqs of length %lu +/- %lu)...", seq_num, seq_len_mean, seq_len_sd);
    ref_genome.generate(seq_num, seq_len_mean, seq_len_sd, ref_nuc_freqs, rng);
  }
  else if (ref_nuc_freqs.size() > 0) {
    // read reference sequence
    fprintf(stderr, "\nReading reference from file '%s'...", fn_ref.c_str());
    ref_genome = Genome(fn_ref.c_str());
  }
  ref_genome.indexRecords();
  fprintf(stderr, "read (%u bp in %u sequences).\n", ref_genome.length, ref_genome.num_records);
  // duplicate genome (all loci homozygous reference)
  ref_genome.duplicate();

  // if a VCF file was provided, apply germline variants
  if (ref_vcf.size() > 0) {
    fprintf(stderr, "applying germline variants (from %s).\n", ref_vcf.c_str());
    VariantSet ref_variants;
    map<string, vector<Genotype >> ref_gt_matrix;
    vario::readVcf(ref_vcf, ref_variants, ref_gt_matrix);
    vario::applyVariants(ref_genome, ref_variants.vec_variants, ref_gt_matrix["0"]);
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
    mutations = vario::generateMutations(num_mutations, random_dbl);
    fprintf(stderr, "\nTotal set of mutations (id, rel_pos, copy):\n");
    for (int i=0; i<num_mutations; i++)
      fprintf(stderr, "%d\t%f\t%d\n", mutations[i].id, mutations[i].relPos, mutations[i].copy);
  }
  // inititalize model of sequence evolution
  SubstitutionModel model;
  if (str_model == "JC") {
    model.init_JC();
  } else if (str_model == "F81") {
    vector<double> vec_freqs = config.getValue<vector<double>>("model-params:nucFreq");
    model.init_F81(&vec_freqs[0]);
  } else if (str_model == "K80") {
    double kappa = config.getValue<double>("model-params:kappa");
    model.init_K80(kappa);
  } else if (str_model == "HKY") {
    vector<double> vec_freqs = config.getValue<vector<double>>("model-params:nucFreq");
    double kappa = config.getValue<double>("model-params:kappa");
    model.init_HKY(&vec_freqs[0], kappa);
  }
  // apply mutations, creating variants in the process
  vector<Variant> variants = vector<Variant>();
  vario::applyMutations(mutations, ref_genome, model, random_dbl, variants);

  // initialize mutation matrix
  vector<vector<short> > mutMatrix(tree.m_numVisibleNodes, std::vector<short>(num_mutations,0));

  // generate clone sequences based on clonal tree and mutations
  //fprintf(stderr, "---\nNow generating clone genomes...\n");
  vector<shared_ptr<Clone>> clones = tree.getVisibleNodes();
  //Clone root = *(tree.m_root);
  for (unsigned i=0; i<clones.size(); ++i) {
    // TODO: this is not how we do things around here, anymore! (update)
    //clones[i]->mutateGenome(ref_genome, mutations, variants, mutMatrix, clone2fn);
  }

  // compile variants for output
  //vector<Clone *> clones = tree.getVisibleNodes();
  int num_nodes = tree.m_numNodes;
  vector<vector<bool> > mat_mut(num_nodes, vector<bool>(num_mutations, false));
  tree.m_root->populateMutationMatrixRec(mat_mut);
  vector<int> vec_vis_nodes_idx = tree.getVisibleNodesIdx();
  vector<string> vec_labels;
  // add node labels
  for (auto i : vec_vis_nodes_idx) {
    vec_labels.push_back(tree.m_vecNodes[i]->label);
  }

  // write clonal variants to file
  ofstream f_vcf;
  f_vcf.open("mutations.vcf");
  vario::writeVcf(ref_genome.records, variants, vec_vis_nodes_idx, vec_labels, mat_mut, f_vcf);
  f_vcf.close();

  // perform ADO
  if (ado_pct > 0)
    for (auto kv : clone2fn)
      seqio::simulateADO(kv.second, ref_genome.masked_length, ado_pct, ado_frag_len, random_dbl);

  return EXIT_SUCCESS;
}
