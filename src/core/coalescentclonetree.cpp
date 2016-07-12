#define DEBUG
#include "coalescentclonetree.hpp"

#include <stdio.h>
#include <vector>


CoalescentCloneTree::CoalescentCloneTree (int numClones) : treeio::Tree<Clone>(numClones) {
#ifdef DEBUG
  _printNodes();
#endif
}

void CoalescentCloneTree::generateRandomTopology(boost::function<float()>& random) {
  // generate N-1 internal nodes (each representing a coalescence event)
  int numClones = m_numNodes;
  int k = numClones-1;
  int nextIndex = numClones;
  for (int i=0; i<numClones-1; i++) {
    // pick first random node (without replacement)
    int index1 = 2*i + random()*k--;
    std::shared_ptr<Clone> p = m_vecNodes[index1];
    m_vecNodes[index1] = m_vecNodes[2*i];
    m_vecNodes[2*i] = p;
#ifdef DEBUG
    fprintf(stderr, "---\nIteration %d:\n", i);
    fprintf(stderr, "\tindex1: %d\n", index1);
    fprintf(stderr, "\tnode1: %s\n", p->label.c_str());
    _printNodes();
#endif
    // pick second random node (without replacement)
    int index2 = 2*i+1 + random()*k;
    std::shared_ptr<Clone> q = m_vecNodes[index2];
    m_vecNodes[index2] = m_vecNodes[2*i+1];
    m_vecNodes[2*i+1] = q;
#ifdef DEBUG
    fprintf(stderr, "\tindex2: %d\n", index2);
    fprintf(stderr, "\tnode2: %s\n", q->label.c_str());
    _printNodes();
#endif
    // create new internal node
    std::shared_ptr<Clone> n(new Clone());
    n->label = ++nextIndex;
    n->m_vecChildren.push_back(p);
    n->m_vecChildren.push_back(q);
    p->parent = n;
    q->parent = n;
    m_vecNodes.push_back(n);
#ifdef DEBUG
    fprintf(stderr, "\tnew internal node: %s\n", n->label.c_str());
    _printNodes();
#endif
  }

  // generate a "healthy" clone as root node
  std::shared_ptr<Clone> r(new Clone());
  r->label = "0";
  r->is_healthy = true;
  r->m_vecChildren.push_back(m_vecNodes[nextIndex-1]);
  r->parent = 0;
  m_vecNodes[nextIndex-1]->parent = r;
  m_vecNodes.push_back(r);
  m_root = r;
#ifdef DEBUG
  fprintf(stderr, "\twe have been ROOTed: %s\n", r->label.c_str());
  _printNodes();
#endif
}

/** Place mutations randomly on the tree.
 * A set of mandatory mutations need to exist between clones,
 * otherwise they cannot be distinguished.
 */
void CoalescentCloneTree::evolve(int numMutations, int numTransformingMutations, boost::function<float()>& random) {
  int num_randomMutations = numMutations - (m_numNodes-1) - numTransformingMutations;
  putMandatoryMutations(numTransformingMutations, random);
  putRandomMutations(num_randomMutations, m_numNodes, random);
}

void CoalescentCloneTree::putMandatoryMutations(int num_transmuts, boost::function<float()>& random) {
fprintf(stderr, "\ndropping %d mandatory mutations (%d transforming + %d to assure distinct clones)...\n", m_numNodes+num_transmuts-1, num_transmuts, m_numNodes-1);
  int m=0;
  // MRCA of clones receives transforming mutations
  std::shared_ptr<Clone> mrca = m_root->m_vecChildren[0];
  while (m<num_transmuts) {
fprintf(stderr, "\tdropping mutation %d on node %s\n", m, mrca->label.c_str());
    mrca->m_vec_mutations.push_back(m);
    m++;
  }
  // each intermediate node needs to introduce a mutation
  _putMandatoryMutation(mrca, m, random);
fprintf(stderr, "%d mutations dropped.\n", m);
}

void CoalescentCloneTree::_putMandatoryMutation(std::shared_ptr<Clone> c, int& mutationId, boost::function<float()>& random) {
  Clone node = *c;
  if (node.isLeaf()) {
    return;
  }
  else {
    bool is_leftMutated = (random()<0.5); // throw a coin
    if (is_leftMutated) {
std::cerr << "\tdropping mutation " << mutationId << " on " << node.m_vecChildren[0] << std::endl;
//fprintf(stderr, "\tdropping mutation %d on node %d\n", mutationId, node.m_vecChildren[0]->label);
      node.m_vecChildren[0]->m_vec_mutations.push_back(mutationId);
    }
    else {
std::cerr << "\tdropping mutation " << mutationId << " on " << node.m_vecChildren[1] << std::endl;
//fprintf(stderr, "\tdropping mutation %d on node %d\n", mutationId, node.m_vecChildren[1]->label);
      node.m_vecChildren[1]->m_vec_mutations.push_back(mutationId);
    }
    mutationId++;
    _putMandatoryMutation(node.m_vecChildren[0], mutationId, random);
    _putMandatoryMutation(node.m_vecChildren[1], mutationId, random);
  }
}

void CoalescentCloneTree::putRandomMutations(int num_mutations, int nextMutationId, boost::function<float()>& random) {
fprintf(stderr, "\nnow dropping %d random mutations...\n", num_mutations);
  while (num_mutations > 0) {
    // mutate random node
    int i = random()*(m_numNodes+1);
    m_vecNodes[i]->m_vec_mutations.push_back(nextMutationId);
fprintf(stderr, "\tdropping mutation %d on node %s\n", nextMutationId, m_vecNodes[i]->label.c_str());
    nextMutationId++;
    num_mutations--;
  }
}

void CoalescentCloneTree::collapseZeroBranches(std::shared_ptr<Clone> node) {
std::cerr << "collapseZeroBranches(" << *node << ")" << std::endl;
//fprintf(stderr, "collapseZeroBranches(<Node(label=%d)>)\n", node->label);
  // check among children whether any has zero distance
  unsigned i=0;
  while ((i<node->m_vecChildren.size()) && ((int)(node->m_vecChildren[i]->distanceToParent())>0)) {
    i++;
  }
  if (i < node->m_vecChildren.size()) { // there is a child with zero distance
    std::shared_ptr<Clone> child = node->m_vecChildren[i];
    // goodbye cruel world!
    Clone::replace(node, child);

std::cerr << *child << " (fresh out of the oven)" << std::endl;
//fprintf(stderr, "<Node(label=%d)> (so fresh & so clean)\n", child->label);
std::cerr << "parent: " << *(child->parent) << std::endl;
//fprintf(stderr, "\tparent <Node(label=%d)>\n", child->parent->label);
for (unsigned k=0; k<child->m_vecChildren.size(); k++) {
std::cerr << "\tchild: " << *(child->m_vecChildren[k]) << std::endl;
//fprintf(stderr, "\tchild <Node(label=%d)>\n", child->m_vecChildren[k]->label);
}
    // repeat process with successor
    if (!child->isLeaf()) {
      collapseZeroBranches(child);
    }
  }
  else { // no immediate child can be collapsed
fprintf(stderr, "No immediate child can be collapsed.\n");
    for (unsigned i=0; i<node->m_vecChildren.size(); i++) {
      collapseZeroBranches(node->m_vecChildren[i]);
    }
  }
}
