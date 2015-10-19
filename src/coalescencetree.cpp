#include "node.hpp"
#include <iostream>
#include <stdio.h>
#include <vector>

/** A binary tree representing the genealogy of samples. */
class CoalescenceTree {
  std::vector<Node> m_vecNodes;
  Node *m_root;

  public:
    CoalescenceTree (int);
    void randomize();
    Node getRoot ();
};

CoalescenceTree::CoalescenceTree (int num_leafs) : m_vecNodes(num_leafs) {
  for (int i=0; i<num_leafs; i++) {
    m_vecNodes[i].label = i;
  }
  m_root = &m_vecNodes[0];
}

Node CoalescenceTree::getRoot () {
  // TODO: return proper root
  return *m_root;
}
