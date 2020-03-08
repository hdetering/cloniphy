#ifndef TREE_H
#define TREE_H

#include "parse.hpp"
#include "../random.hpp"
#include "../clone.hpp"
// #include <functional>

#include <memory>
// #include <cstdio>
#include <string>
#include <vector>

namespace treeio {

namespace parse {
  // forward declaration for later use
  struct node;
}

/** General tree representing a hierarchy of nodes. */
template <typename TNodeType>
struct Tree
{
  std::shared_ptr<TNodeType> m_root;
  int m_numNodes;
  int m_numVisibleNodes;
  int m_numMutations;
  std::vector<std::shared_ptr<TNodeType>> m_vecNodes;

  Tree();
  Tree(int);
  Tree(const parse::node&);
  Tree(const std::string); // read tree from file
  ~Tree();
  /** Return sum of all branch lengths */
  double getTotalBranchLength();
  /** Return vector of relative branch lengths */
  std::vector<double> getAbsoluteBranchLengths();
  /** Return vector of relative branch lengths */
  std::vector<double> getRelativeBranchLengths();
  /** Return visible tree nodes. */
  std::vector<std::shared_ptr<TNodeType>> getVisibleNodes();
  /** Return visible tree nodes' indices. */
  std::vector<int> getVisibleNodesIdx();
  /** Return tree nodes as visited in preorder traversal. */
  std::vector<std::shared_ptr<TNodeType>> getNodesPreOrder();
  void getNodesPreOrderRec(
    const std::shared_ptr<TNodeType> n,
    std::vector<std::shared_ptr<TNodeType>>& nodes);
  /** Return MRCA of non-outgroup lineages */
  std::shared_ptr<TNodeType> getMrca(
    const std::string lbl_outgroup);
  /** Return outgroup lineage */
  std::shared_ptr<TNodeType> getOutgroup(
    const std::string lbl_outgroup);
  /** Build random tree topology. */
  void generateRandomTopology(std::function<double()>&);
  /** Arrange nodes randomly (internal nodes are visible) */
  void generateRandomTopologyInternalNodes(std::function<double()>&);
  /** Arrange nodes randomly (internal nodes are invisible) */
  void generateRandomTopologyLeafsOnly(std::function<double()>&);
  /** Add healthy node as root of clone tree 
   *
   *  This ensures that
   *  a) the healthy genome is represented in the clone tree
   *  b) trunk mutations are assigned to ancestral clone
   */
  void addHealthyRoot (
    const std::string label
  );
  /** Shrink/expand branch length by a random factor
   *  \param lbl_outgrp  Label of outgroup lineage (will not be affected).
   *  \param rnd_dbl     Function returning random double values.
  */
  void 
  varyBranchLengths (
    const std::string lbl_outgrp, 
    std::function<double()>& rand_dbl
  );
  /** Assign random weights to visible nodes */
  void assignWeights(std::vector<double> w);
  /** Distribute somatic mutations along tree branches 
   *  \param n_mut_total   Total number of mutations to drop
   *  \param n_mut_trunk   Number of mutations to place at the most ancestral node
   *  \param lbl_outgroup  Outgroup node label (healthy clone in clone tree)
   */
  virtual void dropSomaticMutations(
    const int n_mut_total,
    const int n_mut_trunk,
    const std::string lbl_outgroup,
    RandomNumberGenerator& rng);
  void printNewick(const std::string);
  void printNewick(std::ostream&);
  void printNewick(std::shared_ptr<TNodeType>, std::ostream&, bool first=true);
  void printNexus(const std::string);
  void printNexus(std::ostream&);
  void printNexusTranslate(std::ostream&);
  void printNexusMatrix(std::ostream&);
  void printNexusTree(std::shared_ptr<TNodeType>, std::ostream&, bool first=true);
  void printDot(const std::string filename);
  void printDot(std::shared_ptr<TNodeType>, std::ostream&);
  /** outputs boolean matrix of mutational states for visible nodes */
  void writeMutationMatrix(const std::string filename);
  /** outputs boolean matrix of mutational states for visible nodes */
  void writeMutationMatrix(std::ostream& outstream);

  // useful for debugging purposes:
  void _printNodes();
  void _printTreeInfo();

private:
  /** Assign initial mutations to founding clone.
   *  \param sp_mrca  MRCA of non-outgroup lineages.
   *  \param n_muts   Number of mutations to assign.
   */
  void dropTransformingMutations(
    std::shared_ptr<TNodeType> sp_mrca,
    const int n_muts
  );
  /** Make sure each clone has at least 1 mutation difference to every other clone. */
  void dropMandatoryMutations(std::shared_ptr<TNodeType>, int&);
  /** Drop free mutations on random clone nodes. 
   *  \param lbl_outgroup  Label of outgroup lineage.
   *  \param n_muts        Number of mutations to assign.
   *  \param id_mut        Identifier of next mutation (for consistent numbering).
   *  \param rng           Random number generator.
   */
  void dropRandomMutations(
    const std::string lbl_outgroup,
    const int n_muts,
    int& id_mut, 
    RandomNumberGenerator& rng);
  /** Reset mutation ids to follow pre-oder traversal sequence. */
  void _relabelMutationsRec(std::shared_ptr<TNodeType>, int&);
  void _varyBranchLengthsRec(std::shared_ptr<TNodeType>, std::function<double()>&);
  void _assignWeightsRec(
    std::shared_ptr<TNodeType>,
    std::vector<double>::iterator&,
    std::vector<double>);
  void _printDotRec(std::shared_ptr<TNodeType>, std::ostream&);
  std::shared_ptr<TNodeType> _adaptNode(const parse::node&);
};

} // namespace treeio

#endif // TREE_H