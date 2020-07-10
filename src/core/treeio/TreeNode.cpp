#include "TreeNode.hpp"
#include <iostream>

using namespace std;

namespace treeio {

TreeNode::TreeNode() {
  parent = 0;
}

TreeNode::~TreeNode() {}

// streaming operator for easy printing
ostream& operator<<(std::ostream& stream, const TreeNode& node) {
  stream << "<TreeNode(index=" << node.index << ", label='" << node.label << "', length=" << node.length << ")>";
  return stream;
}

bool TreeNode::isRoot() {
  if (this->parent == 0)
    return true;
  else
    return false;
}

} // namespace treeio