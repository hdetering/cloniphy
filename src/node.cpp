#include "node.hpp"
#include <cstddef>

Node::Node() {
  parent = NULL;
  left = NULL;
  right = NULL;
}

bool Node::isLeaf() {
  return ((left == NULL) && (right == NULL));
}

bool Node::isRoot() {
  return (parent == NULL);
}
