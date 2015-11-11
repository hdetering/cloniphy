#ifndef COALESCENTCLONETREE_H
#define COALESCENTCLONETREE_H

#include "clonetree.hpp"
#include "clone.hpp"
#include <boost/function.hpp>
#include <vector>

/** Binary tree representing the genealogy of samples. */
class CoalescentCloneTree : public CloneTree
{
  protected:
    /** Make sure each clone has at least 1 mutation difference to every other clone. */
    void putMandatoryMutations(int, boost::function<float()>&);
    /** Drop mutations recursively. */
    void _putMandatoryMutation(Clone *c, int& mutationId, boost::function<float()>& random);
    /** Drop free mutations on random clone nodes. */
    void putRandomMutations(int, int, boost::function<float()>&);

  public:
    CoalescentCloneTree(int, std::vector<float>);
    /** Build random coalescence tree containing clones as tips. */
    void generateRandomTopology(boost::function<float()>&);
    /** Drop random mutations on clones. */
    void evolve(int, int, boost::function<float()>&);
    /** Replace clones that have identical children. */
    void collapseZeroBranches(Clone *);
};

#endif /* COALESCENTCLONETREE_H */
