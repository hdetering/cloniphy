#ifndef CLONETREE_H
#define CLONETREE_H

#include "clone.hpp"

/** General tree representing a hierarchy of clones. */
class CloneTree {
  protected:
    Clone *m_root;
    int m_numClones;
    std::vector<Clone *> m_vecNodes;

  public:
    CloneTree(int, std::vector<float>);
    /** Returns the clone tree's root node. */
    Clone* getRoot();
};

#endif /* CLONETREE_H */
