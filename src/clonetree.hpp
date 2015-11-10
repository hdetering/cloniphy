#ifndef CLONETREE_H
#define CLONETREE_H

#include "clone.hpp"

/** General tree representing a hierarchy of clones. */
class CloneTree {
  protected:
    Clone *m_root;
    int m_numClones;
    std::vector<Clone *> m_vecNodes;

    /** Recursive part of generating DOT representation. */
    static void _printDotRecursive(Clone *, std::ostream&);

  public:
    CloneTree(int, std::vector<float>);
    /** Returns the clone tree's root node. */
    Clone* getRoot();
    /** Returns all visible nodes of the clone tree. */
    std::vector<Clone *> getVisibleNodes();
    /** Generates the graph representation of a (sub)tree in DOT format. */
    static void printDot(Clone *, std::ostream& = std::cout);
};

#endif /* CLONETREE_H */
