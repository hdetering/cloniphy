#include "clonetree.hpp"

/** C'tor creates a vector of clones and assignes frequency to each clone. */
CloneTree::CloneTree (int numClones, std::vector<float> freqs) : m_numClones(numClones), m_vecNodes(numClones) {
  // initialize tips
  for (int i=0; i<numClones; i++) {
    Clone *c = new Clone();
    c->label = i+1;
    c->freq = freqs[i];
    c->is_visible = true;
    m_vecNodes[i] = c;
  }
#ifdef DEBUG
  printNodes();
#endif
}

Clone* CloneTree::getRoot() {
  return m_root;
}
