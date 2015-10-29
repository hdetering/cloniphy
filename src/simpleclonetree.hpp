#ifndef SIMPLECLONETREE_H
#define SIMPLECLONETREE_H

#include "tree.hpp"
#include "clone.hpp"
#include <boost/function.hpp>
#include <vector>

/** Binary tree representing the genealogy of samples. */
class SimpleCloneTree : public Tree
{
  protected:
    int m_numClones;
    Clone *m_root;
    std::vector<Clone *> m_vecNodes;
    std::vector<Clone *> m_vecLeafs;

    /** Make sure each clone has at least 1 mutation difference to every other clone. */
    void putMandatoryMutations(int, boost::function<float()>&);
    /** Drop mutations recursively. */
    void _putMandatoryMutation(Clone *c, int& mutationId, boost::function<float()>& random);
    /** Drop free mutations on random clone nodes. */
    void putRandomMutations(int, int, boost::function<float()>&);
    /** Recursive part of generating DOT representation. */
    static void _printDotRecursive(Clone *, std::ostream&);
    /** Recursive part of generating Newick representation. */
    static void _printNewickRecursive(Clone *, bool, std::ostream&);

  public:
    SimpleCloneTree(int, std::vector<float>);
    /** Returns the clone tree's root node. */
    Clone* getRoot();
    /** Returns visible clones. */
    std::vector<Clone *> getLeafs();
    /** Build random coalescence tree containing clones as tips. */
    void generateRandomTopology(boost::function<float()>&);
    /** Drop random mutations on clones. */
    void evolve(int, boost::function<float()>&);
    /** Replace clones that have identical children. */
    void collapseZeroBranches(Clone *);
    /** Print vector containing clones to STDERR (for debugging). */
    void printNodes();
    /** Generates the graph representation of a (sub)tree in DOT format. */
    static void printDot(Clone *, std::ostream& = std::cout);
    /** Generates the string representation of a (sub)tree in Newick notation. */
    static void printNewick(Clone *, std::ostream& = std::cout);
};

#endif /* SIMPLECLONETREE_H */
