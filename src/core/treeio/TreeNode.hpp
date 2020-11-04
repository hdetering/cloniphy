#ifndef TREENODE_H
#define TREENODE_H

// #include "random.hpp"
// #include <functional>

#include <memory>
// #include <cstdio>
#include <string>
#include <vector>

namespace treeio {

/** Generic tree node */
struct TreeNode
{
  int index;
  std::string label;
  double length;
  double weight;
  bool is_visible;
  std::shared_ptr<TreeNode> parent;
  std::vector<std::shared_ptr<TreeNode>> m_vecChildren;

  TreeNode();
  virtual ~TreeNode();
  virtual float distanceToParent() = 0;
  virtual bool isLeaf() = 0;
  bool isRoot();
};

// streaming operator for easy printing
std::ostream& operator<<(std::ostream&, const TreeNode&);

} // namespace treeio

#endif // TREENODE_H