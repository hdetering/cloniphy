#ifndef CLONE_H
#define CLONE_H

#include "seqio.hpp"
#include "treeio.hpp"
#include "vario.hpp"
#include "evolution.hpp"
#include <functional>
#include <map>
#include <memory>
#include <ostream>
#include <vector>

using seqio::Genome;
using vario::Genotype;
using vario::Mutation;
using vario::Variant;
using evolution::GermlineSubstitutionModel;

struct Clone: public treeio::TreeNode
{
  float freq;
  bool is_healthy;
  std::vector<std::shared_ptr<Clone>> m_vecChildren;
  std::vector<int> m_vec_mutations;
  std::vector<Genotype> m_vec_genotypes;

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
