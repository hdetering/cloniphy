#ifndef SIMPLECLONETREE_H
#define SIMPLECLONETREE_H

#include "tree.hpp"
#include "clone.hpp"
#include <vector>

/** Binary tree representing the genealogy of samples. */
class SimpleCloneTree : public Tree
{
  public:
    SimpleCloneTree(int, std::vector<float>);
    void generateRandomTopology(long);
};

#endif /* SIMPLECLONETREE_H */
