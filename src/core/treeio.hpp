#ifndef TREEIO_H
#define TREEIO_H

#include <functional>
#include <iostream>
#include <memory>
#include <stdio.h>
#include <string>
#include <vector>

namespace treeio {

/** Generic tree node */
struct TreeNode
{
  int index;
  std::string label;
  double length;
  double weight;
  bool is_visible;
  std::shared_ptr<TreeNode> parent;
  std::vector<std::shared_ptr<TreeNode>> m_vecChildren;

  TreeNode();
  virtual ~TreeNode();
  virtual float distanceToParent() = 0;
  virtual bool isLeaf() = 0;
  bool isRoot();
};

// streaming operator for easy printing
std::ostream& operator<<(std::ostream&, const TreeNode&);

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
  /** Build random tree topology. */
  void generateRandomTopology(std::function<double()>&);
  /** Arrange nodes randomly (internal nodes are visible) */
  void generateRandomTopologyInternalNodes(std::function<double()>&);
  /** Arrange nodes randomly (internal nodes are invisible) */
  void generateRandomTopologyLeafsOnly(std::function<double()>&);
  /** Shrink/expand branch length by a random factor */
  void varyBranchLengths(std::function<double()>&);
  /** Assign random weights to visible nodes */
  void assignWeights(std::vector<double> w);
  /** Distribute somatic mutations along tree branches */
  virtual void dropSomaticMutations(
    int n_mut_total,
    int n_mut_trunk,
    RandomNumberGenerator<>& rng);
  void printNewick(const std::string);
  void printNewick(std::ostream&);
  void printNewick(std::shared_ptr<TNodeType>, std::ostream&, bool first=true);
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
  /** Assign initial mutations to founding clone. */
  void dropTransformingMutations(int);
  /** Make sure each clone has at least 1 mutation difference to every other clone. */
  void dropMandatoryMutations(std::shared_ptr<TNodeType>, int&);
  /** Drop free mutations on random clone nodes. */
  void dropRandomMutations(int, int&, RandomNumberGenerator<>&);
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

/** Reads and writes tree files. */
namespace parse {

// typedef to ease the writing
typedef std::vector<node> children_vector;

struct node
{
  std::string label;
  double length = 0;
  children_vector children;
};

void readNewick(std::string, node&);
void readNewick(std::istream&, node&);

} // namespace parse

} // namespace treeio

#endif /* TREEIO_H */
