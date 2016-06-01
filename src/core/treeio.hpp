#ifndef TREEIO_H
#define TREEIO_H

#include <boost/function.hpp>
#include <iostream>
#include <stdio.h>
#include <string>
#include <vector>

/** Generic tree node */
struct Node
{
  int index;
  std::string label;
  double length;
  bool is_visible;
  Node *parent;
  std::vector<Node *> m_vecChildren;

  Node();
  virtual ~Node();
  virtual float distanceToParent() = 0;
  virtual bool isLeaf() = 0;
  bool isRoot();
};

// streaming operator for easy printing
std::ostream& operator<<(std::ostream&, const Node&);

namespace treeio {
  // forward declaration for later use
  struct node;
}

/** General tree representing a hierarchy of nodes. */
template <typename T>
struct Tree
{
  T *m_root;
  int m_numNodes;
  std::vector<T *> m_vecNodes;

  Tree();
  Tree(int);
  Tree(const treeio::node&);
  ~Tree();
  /** Return sum of all branch lengths */
  double getTotalBranchLength();
  /** Return vector of relative branch lengths */
  std::vector<double> getRelativeBranchLengths();
  /** Return visible tree nodes. */
  std::vector<T *> getVisibleNodes();
  /** Build random tree topology. */
  void generateRandomTopology(boost::function<double()>&);
  /** Arrange nodes randomly (internal nodes are visible) */
  void generateRandomTopologyInternalNodes(boost::function<double()>&);
  /** Arrange nodes randomly (internal nodes are invisible) */
  void generateRandomTopologyLeafsOnly(boost::function<double()>&);
  /** Drop random mutations on clones. */
  virtual void evolve(int, int, boost::function<double()>&);
  void printDot(T*, std::ostream&);

private:
  /** Assign initial mutations to founding clone. */
  void dropTransformingMutations(int);
  /** Make sure each clone has at least 1 mutation difference to every other clone. */
  void dropMandatoryMutations(T*, int&);
  /** Drop free mutations on random clone nodes. */
  void dropRandomMutations(int, int&, boost::function<double()>&);
  void _printDotRecursive(T*, std::ostream&);
  void _printNodes();
  T* _adaptNode(const treeio::node&);
};


/** Reads and writes tree files. */
namespace treeio {
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
}

#endif /* TREEIO_H */
