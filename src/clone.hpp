#ifndef CLONE_H
#define CLONE_H

#include "seqio.hpp"
#include "treeio.hpp"
#include "vario.hpp"
#include "evolution.hpp"
#include <boost/function.hpp>
#include <ostream>
#include <vector>

using seqio::Genome;
using vario::Genotype;
using vario::Mutation;
using vario::Variant;
using evolution::SubstitutionModel;

struct Clone: public Node
{
  float freq;
  bool is_healthy;
  Clone *parent;
  std::vector<Clone *> m_vecChildren;
  std::vector<int> m_vecMutations;
  std::vector<Variant> m_vec_variants;
  std::vector<Genotype> m_vec_genotypes;

  Clone();
  virtual ~Clone();
  void setParent(Clone *);
  bool isLeaf();
  std::vector<Clone *> getChildren();
  float distanceToParent();
  /** replace other clone in the tree (needed to collapse branches) */
  void replace(Clone *);
  /** modify the given sequence by applying a set of mutations. */
  void mutateGenome(Genome&, const std::vector<Mutation>&,
                    std::vector<Variant>&, SubstitutionModel,
                    std::vector<std::vector<short> >&,
                    boost::function<float()>& random);
};

#endif /* CLONE_H */
