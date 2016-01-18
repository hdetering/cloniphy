#ifndef BASICCLONETREE_H
#define BASICCLONETREE_H

#include "clonetree.hpp"
#include "clone.hpp"
#include <boost/function.hpp>
#include <vector>

class BasicCloneTree : public CloneTree {
  protected:
    /** Drop cancer-initiating mutations on top clones. */
    void dropTransformingMutations(int);
    /** Drop one mutation on a given node and repeat for children. */
    void dropMandatoryMutations(Clone*, int&);
    /** Drop mutations randomly onto clone nodes. */
    void dropRandomMutations(int, int&, boost::function<float()>&);

  public:
    BasicCloneTree();
    BasicCloneTree(int);
    BasicCloneTree(const treeio::node&);

    /** Build random clone tree. */
    void generateRandomTopology(boost::function<float()>&);
    /** Drop random mutations on clones. */
    void evolve(int, int, boost::function<float()>&);
};

#endif /* BASICCLONETREE_H */
