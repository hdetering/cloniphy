#ifndef CLONE_H
#define CLONE_H

#include "treeio/TreeNode.hpp"
#include <functional>
#include <map>
#include <memory>
#include <ostream>
#include <vector>

struct Clone: public treeio::TreeNode
{
  bool is_healthy;
  std::vector<std::shared_ptr<Clone>> m_vecChildren;
  std::vector<int> m_vec_mutations;

  Clone();
  virtual ~Clone();
  static void setParent(std::shared_ptr<Clone> child, std::shared_ptr<Clone> parent);
  bool isLeaf();
  std::vector<std::shared_ptr<Clone>> getChildren();
  /** Fill a matrix (nodes x mutations) denoting mutational state of subtree rooted in this node. */
  void populateMutationMatrixRec(std::vector<std::vector<bool>> &m);
  /** Fill a matrix (nodes x mutations) denoting nodes in which mutations first occurred. */
  void populateMutationOccRec(std::vector<std::vector<bool>> &m);
  float distanceToParent();
  /** replace other clone in the tree (needed to collapse branches) */
  static void replace(std::shared_ptr<Clone> old_node, std::shared_ptr<Clone> new_node);
};

#endif /* CLONE_H */
