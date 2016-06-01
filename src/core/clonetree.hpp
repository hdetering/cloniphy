#ifndef CLONETREE_H
#define CLONETREE_H

#include "clone.hpp"
#include "treeio.hpp"
#include <boost/function.hpp>
#include <iostream>
#include <vector>

/** General tree representing a hierarchy of clones. */
class CloneTree {
  protected:
    Clone *m_root;
    int m_numClones;
    std::vector<Clone *> m_vecNodes;

    /** Recursive part of generating Newick representation. */
    static void _printNewickRecursive(Clone *, bool, std::ostream&);
    /** Recursive part of generating DOT representation. */
    static void _printDotRecursive(Clone *, std::ostream&);

  public:
    CloneTree();
    CloneTree(int);
    CloneTree(const treeio::node&);
    ~CloneTree();

    /* --------------------------------------*
     *          Virtual functions            *
     * Need to be implemented by subclasses. *
     *---------------------------------------*/

    /** Build random tree topology. */
    virtual void generateRandomTopology(boost::function<float()>&) = 0;
    /** Drop random mutations on clones. */
    virtual void evolve(int, int, boost::function<float()>&) = 0;

    /* ----------------------------------*
     * Functions shared by clonal trees. *
     * ----------------------------------*/

    /** Returns the clone tree's root node. */
    Clone* getRoot();
    /** Returns all visible nodes of the clone tree. */
    std::vector<Clone *> getVisibleNodes();
    /** Generates the string representation of a (sub)tree in Newick notation. */
    static void printNewick(Clone *, std::ostream& = std::cout);
    /** Generates the graph representation of a (sub)tree in DOT format. */
    static void printDot(Clone *, std::ostream& = std::cout);
    /** Print vector containing clones to STDERR (for debugging). */
    void printNodes();
    /** Initialize CloneTree from generic tree */
    Clone* adaptFromGeneric(const treeio::node);
};

#endif /* CLONETREE_H */
