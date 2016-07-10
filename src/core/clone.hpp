#ifndef CLONE_H
#define CLONE_H

#include "seqio.hpp"
#include "treeio.hpp"
#include "vario.hpp"
#include "evolution.hpp"
#include <functional>
#include <map>
#include <ostream>
#include <vector>

using seqio::Genome;
using vario::Genotype;
using vario::Mutation;
using vario::Variant;
using evolution::SubstitutionModel;

struct Clone: public treeio::TreeNode
{
  float freq;
  bool is_healthy;
  Clone *parent;
  std::vector<Clone *> m_vecChildren;
  std::vector<int> m_vec_mutations;
  std::vector<Genotype> m_vec_genotypes;

  Clone();
  virtual ~Clone();
  void setParent(Clone *);
  bool isLeaf();
  std::vector<Clone *> getChildren();
  /** Fill a matrix (nodes x mutations) denoting mutational state of subtree rooted in this node. */
  void populateMutationMatrixRec(std::vector<std::vector<bool>> &m);
  /** Fill a matrix (nodes x mutations) denoting nodes in which mutations first occurred. */
  void populateMutationOccRec(std::vector<std::vector<bool>> &m);
  float distanceToParent();
  /** replace other clone in the tree (needed to collapse branches) */
  void replace(Clone *);
  /** create a modified genome based on given sequence by applying a set of mutations. */
  void mutateGenome(
    const Genome&,
    const std::vector<Mutation>&,
    std::vector<Variant>&,
    std::vector<std::vector<short> >&,
    std::map<Clone*, std::string>&);
  /** apply a predefined set of mutations to own genome */
  void applyMutations(
    const std::vector<Mutation>&,
    SubstitutionModel,
    Genome&,
    std::function<double()>&);
};

#endif /* CLONE_H */
