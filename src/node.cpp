#include "node.hpp"
#include <cstddef>

Node::Node() {
  parent = NULL;
}

Node::~Node() {}

bool Node::isRoot() {
  return (parent == NULL);
}
