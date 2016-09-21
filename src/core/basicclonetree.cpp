#include "basicclonetree.hpp"
#include <stdio.h>
#include <vector>

BasicCloneTree::BasicCloneTree() : treeio::Tree<Clone>() {}
BasicCloneTree::BasicCloneTree(int numClones) : treeio::Tree<Clone>(numClones) {}
BasicCloneTree::BasicCloneTree(const treeio::parse::node& root) : treeio::Tree<Clone>(root) {}

void BasicCloneTree::generateRandomTopology(boost::function<float()>& random) {
  // generate a "healthy" clone as root node
  std::shared_ptr<Clone> r(new Clone());
  r->label = "0";
  r->is_healthy = true;
  r->parent = 0;
  m_vecNodes.push_back(r);
  m_root = r;
  // first clone becomes child of root
  Clone::setParent(m_vecNodes[0], r);

  // pick parent for each clone
  std::vector<std::shared_ptr<Clone>> parents;
  parents.push_back(m_vecNodes[0]);
  for (int i=1; i<m_numNodes; ++i) {
    std::shared_ptr<Clone> c = m_vecNodes[i];
    int p_index = random()*i;
    std::shared_ptr<Clone> p = parents[p_index];
std::cerr << *c << " gets parent " << *p << std::endl;
//fprintf(stderr, "\tClone<label=%s> gets parent Clone<label=%s>\n", c->label.c_str(), p->label.c_str());
    Clone::setParent(c, p);
    parents.push_back(c);
  }
}

void BasicCloneTree::evolve(int numMutations, int numTransformingMutations, boost::function<float()>& random) {
fprintf(stderr, "Dropping %d mutations (%d transforming)...\n", numMutations, numTransformingMutations);
  dropTransformingMutations(numTransformingMutations);
  int nextMutId = numTransformingMutations; // identifier for next mutation
fprintf(stderr, "Now dropping mandatory mutations...\n");
  // each clone needs at least one private mutation
  dropMandatoryMutations(m_root, nextMutId);
  // remaining mutations are dropped randomly
  int numRandomMutations = numMutations - numTransformingMutations - m_numNodes;
fprintf(stderr, "Now dropping %d random mutations...\n", numRandomMutations);
  dropRandomMutations(numRandomMutations, ++nextMutId, random);
}

/** Drop transforming mutations on immediate children of root node. */
void BasicCloneTree::dropTransformingMutations(int numMutations) {
  std::vector<std::shared_ptr<Clone>> topNodes = m_root->m_vecChildren;
  for (unsigned i=0; i<topNodes.size(); ++i) {
    std::shared_ptr<Clone> c = topNodes[i];
    for (int m=0; m<numMutations; ++m) {
std::cerr << "\tDropping mutation " << m << " on " << *c << std::endl;;
//fprintf(stderr, "\tDropping mutation %d on Clone<label=%s>\n", m, c->label.c_str());
      c->m_vec_mutations.push_back(m);
    }
  }
}

/** Drop one mutation on a given node and repeat for children. */
void BasicCloneTree::dropMandatoryMutations(std::shared_ptr<Clone> clone, int &mutationId) {
  if (clone!=m_root) {
std::cerr << "\tDropping mutation " << mutationId << " on " << *clone << std::endl;
//fprintf(stderr, "\tDropping mutation %d on Clone<label=%d>\n", mutationId, clone->label);
    clone->m_vec_mutations.push_back(mutationId);
  }
  else {
    mutationId--;
  }
  for (unsigned i=0; i<clone->m_vecChildren.size(); ++i) {
    dropMandatoryMutations(clone->m_vecChildren[i], ++mutationId);
  }
}

void BasicCloneTree::dropRandomMutations(int numMutations, int &mutationId, boost::function<float()>& random) {
  for (int i=0; i<numMutations; ++i) {
    // pick random clone to mutate
    std::shared_ptr<Clone> c = m_vecNodes[random()*m_numNodes];
std::cerr << "\tDropping mutation " << mutationId << " on " << *c << std::endl;
//fprintf(stderr, "\tDropping mutation %d on Clone<label=%d>\n", mutationId, c->label);
    c->m_vec_mutations.push_back(mutationId++);
  }
}
