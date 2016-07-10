#include "clonetree.hpp"
#include "treeio.hpp"
#include <boost/lexical_cast.hpp>
#include <stdio.h>

CloneTree::CloneTree() {}

/** C'tor creates a vector of clones and assignes frequency to each clone. */
CloneTree::CloneTree(int numClones) : m_numClones(numClones), m_vecNodes(numClones) {
  // initialize tips
  for (int i=0; i<numClones; i++) {
    Clone *c = new Clone();
    c->index = i+1;
    c->label = boost::lexical_cast<std::string>(i+1);
    //c->freq = freqs[i];
    c->is_visible = true;
    m_vecNodes[i] = c;
  }
#ifdef DEBUG
  printNodes();
#endif
}
CloneTree::CloneTree(const treeio::parse::node& root) {
  this->m_numClones = 0;
  Clone *root_clone = adaptFromGeneric(root);
  this->m_root = root_clone;
}
CloneTree::~CloneTree() {
  for (unsigned i=0; i<m_vecNodes.size(); ++i) {
     delete(m_vecNodes[i]);
  }
}

Clone* CloneTree::adaptFromGeneric(const treeio::parse::node node) {
  Clone *clone = new Clone();
  clone->index = m_numClones++;
  clone->label = node.label;
  clone->length = node.length;
  clone->is_visible = true;
  this->m_vecNodes.push_back(clone);

  for (unsigned i=0; i<node.children.size(); ++i) {
    Clone *child_clone = adaptFromGeneric(node.children[i]);
    clone->m_vecChildren.push_back(child_clone);
    child_clone->parent = clone;
  }

  return clone;
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

/** Generates the string representation of a (sub)tree in Newick notation. */
void CloneTree::printNewick(Clone *node, std::ostream& os) {
  _printNewickRecursive(node, true, os);
  os << ";";
}

/** Prints the string representation of the tree in Newick notation. (internal) */
void CloneTree::_printNewickRecursive(Clone *node, bool isFirstChild, std::ostream& os) {
  if (!isFirstChild) { os << ","; }
  if (node->getChildren().size() > 0) {
    os << "(";
    isFirstChild = true;
    for (unsigned i=0; i<node->getChildren().size(); i++) {
      _printNewickRecursive(node->getChildren()[i], isFirstChild, os);
      isFirstChild = false;
    }
    os << ")" << node->label << ":" << node->distanceToParent()+0.05;
  }
  else {
    os << node->label << ":" << node->distanceToParent()+0.05;
  }
}

void CloneTree::printDot(Clone *node, std::ostream& os) {
  os << "digraph G {\n";
  _printDotRecursive(node, os);
  os << "}\n";
}

void CloneTree::_printDotRecursive(Clone *node, std::ostream& os) {
  if (node->isRoot()) {
    os << "\t" << node->index << " [style=filled,color=limegreen];" << std::endl;
  } else if (node->is_visible) {
    os << "\t" << node->index << " [style=filled,color=tomato];" << std::endl;
  }
  for (unsigned i=0; i<node->m_vecChildren.size(); ++i) {
    Clone *child = node->m_vecChildren[i];
    float edgeWeight = child->distanceToParent();
    os << "\t" << node->index << " -> " << child->index;
    if (edgeWeight > 0.0) {
      os << "[style=bold,label=" << edgeWeight << "]";
    }
    os << ";" << std::endl;
    os << "\t" << node->index << " [label=\"" << node->label << "\"];" << std::endl;

    _printDotRecursive(child, os);
  }
}

// use me for debugging :-)
void CloneTree::printNodes() {
  for (unsigned i=0; i<m_vecNodes.size(); i++) { fprintf(stderr, "|%2u ", i); }; fprintf(stderr, "|\n");
  for (unsigned i=0; i<m_vecNodes.size(); i++) { fprintf(stderr, "+---"); }; fprintf(stderr, "+\n");
  for (unsigned i=0; i<m_vecNodes.size(); i++) {
    if (m_vecNodes[i]) { fprintf(stderr, "|%s ", m_vecNodes[i]->label.c_str()); }
    else { fprintf(stderr, "| - "); }
  };
  fprintf(stderr, "|\n");
}
