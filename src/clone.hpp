#ifndef CLONE_H
#define CLONE_H

#include "node.hpp"
#include "seqio.hpp"
#include <ostream>
#include <vector>

struct Clone: public Node
{
  float freq;
  bool is_healthy;
  bool is_visible;
  Clone *parent;
  std::vector<Clone*> m_vecChildren;
  std::vector<int> m_vecMutations;

  Clone();
  virtual ~Clone();
  bool isLeaf();
  std::vector<Clone *> getChildren();
  float distanceToParent();
  /** replace other clone in the tree (needed to collapse branches) */
  void replace(Clone *);
  /** modify the given sequence by applying a set of mutations. */
  void mutateGenome(std::vector<SeqRecord>&, const std::vector<long>&, const std::vector<Mutation>&);
};

#endif /* CLONE_H */
