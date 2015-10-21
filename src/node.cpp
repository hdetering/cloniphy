#include "node.hpp"
#include <cstddef>

Node::Node() {
  parent = NULL;
}

bool Node::isLeaf() {
  return m_vecChildren.empty();
}

bool Node::isRoot() {
  return (parent == NULL);
}
