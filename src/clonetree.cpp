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

std::vector<Clone *> CloneTree::getVisibleNodes() {
  std::vector<Clone *> vis_nodes;
  for (unsigned i=0; i<m_vecNodes.size(); ++i) {
    Clone *c = m_vecNodes[i];
    if (c->is_visible && (c!=m_root)) {
      vis_nodes.push_back(c);
    }
  }
  return vis_nodes;
}

void CloneTree::printDot(Clone *node, std::ostream& os) {
  os << "digraph G {\n";
  _printDotRecursive(node, os);
  os << "}\n";
}

void CloneTree::_printDotRecursive(Clone *node, std::ostream& os) {
  if (node->is_healthy) {
    os << "\t" << node->label << " [style=filled,color=limegreen];" << std::endl;
  } else if (node->is_visible) {
    os << "\t" << node->label << " [style=filled,color=tomato];" << std::endl;
  }
  for (unsigned i=0; i<node->m_vecChildren.size(); ++i) {
    Clone *child = node->m_vecChildren[i];
    float edgeWeight = child->distanceToParent();
    os << "\t" << node->label << " -> " << child->label;
    if (edgeWeight > 0.0) {
      os << "[style=bold,label=" << edgeWeight << "]";
    }
    os << ";" << std::endl;

    _printDotRecursive(child, os);
  }
}
