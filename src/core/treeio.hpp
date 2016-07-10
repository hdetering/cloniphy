#ifndef TREEIO_H
#define TREEIO_H

#include <functional>
#include <iostream>
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
  TreeNode *parent;
  std::vector<TreeNode *> m_vecChildren;

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
template <typename T>
struct Tree
{
  T *m_root;
  int m_numNodes;
  int m_numVisibleNodes;
  std::vector<T *> m_vecNodes;

  Tree();
  Tree(int);
  Tree(const parse::node&);
  ~Tree();
  /** Return sum of all branch lengths */
  double getTotalBranchLength();
  /** Return vector of relative branch lengths */
  std::vector<double> getAbsoluteBranchLengths();
  /** Return vector of relative branch lengths */
  std::vector<double> getRelativeBranchLengths();
  /** Return visible tree nodes. */
  std::vector<T *> getVisibleNodes();
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
  /** Drop random mutations on clones. */
  virtual void evolve(int, int, RandomNumberGenerator<>&);
  void printNewick(std::ostream&);
  void printNewick(T*, std::ostream&, bool first=true);
  void printDot(T*, std::ostream&);

private:
  /** Assign initial mutations to founding clone. */
  void dropTransformingMutations(int);
  /** Make sure each clone has at least 1 mutation difference to every other clone. */
  void dropMandatoryMutations(T*, int&);
  /** Drop free mutations on random clone nodes. */
  void dropRandomMutations(int, int&, RandomNumberGenerator<>&);
  void _varyBranchLengthsRec(T*, std::function<double()>&);
  void _assignWeightsRec(T*, std::vector<double>::iterator&, std::vector<double>);
  void _printDotRec(T*, std::ostream&);
  void _printNodes();
  T* _adaptNode(const parse::node&);
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
