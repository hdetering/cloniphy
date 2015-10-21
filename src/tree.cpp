#include "tree.hpp"
#include <iostream>

Tree::Tree(int numTips) : m_numTips(numTips), m_vecNodes(2*numTips-1) {
  // initialize pointers
  for (int i=0; i<numTips; i++) {
    m_vecNodes[i] = 0;
  }
}

Node* Tree::getRoot() {
  return m_root;
}

// use me for debugging :-)
void Tree::printNodes() {
  for (unsigned i=0; i<m_vecNodes.size(); i++) { fprintf(stderr, "|%2u ", i); }; fprintf(stderr, "|\n");
  for (unsigned i=0; i<m_vecNodes.size(); i++) { fprintf(stderr, "+---"); }; fprintf(stderr, "+\n");
  for (unsigned i=0; i<m_vecNodes.size(); i++) {
    if (m_vecNodes[i] > 0) { fprintf(stderr, "|%2d ", m_vecNodes[i]->label); }
    else { fprintf(stderr, "| - "); }
  };
  fprintf(stderr, "|\n");
}

/** Prints the string representation of the tree in Newick notation. */
void Tree::printNewick(Node *node, std::ostream& os) {
  _printNewick(node, true, os);
  os << ";";
}

/** Prints the string representation of the tree in Newick notation. (internal) */
void Tree::_printNewick(Node *node, bool isFirstChild, std::ostream& os) {
  Node n = *node;

  if (!isFirstChild) { os << ","; }
  if (n.m_vecChildren.size() > 0) {
    os << "(";
    isFirstChild = true;
    for (unsigned i=0; i<n.m_vecChildren.size(); i++) {
      _printNewick(n.m_vecChildren[i], isFirstChild, os);
      isFirstChild = false;
    }
    os << ")" << n.label;
  }
  else {
    os << n.label;
  }
}
