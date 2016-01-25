#ifndef CLONE_H
#define CLONE_H

#include "seqio.hpp"
#include "treeio.hpp"
#include "vario.hpp"
#include <ostream>
#include <vector>

using seqio::SeqRecord;
using vario::Mutation;

struct Clone: public Node
{
  float freq;
  bool is_healthy;
  Clone *parent;
  std::vector<Clone *> m_vecChildren;
  std::vector<int> m_vecMutations;

  Clone();
  virtual ~Clone();
  void setParent(Clone *);
  bool isLeaf();
  std::vector<Clone *> getChildren();
  float distanceToParent();
  /** replace other clone in the tree (needed to collapse branches) */
  void replace(Clone *);
  /** modify the given sequence by applying a set of mutations. */
  void mutateGenome(std::vector<SeqRecord>&, const std::vector<unsigned long>&, const std::vector<Mutation>&, std::vector<short>&);
};

#endif /* CLONE_H */
