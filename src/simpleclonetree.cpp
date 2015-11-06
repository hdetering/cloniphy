#define DEBUG
#include "simpleclonetree.hpp"
#include <boost/random.hpp>
// Choosing the random number generator. (mt19937: Mersenne-Twister)
typedef boost::mt19937 base_generator_type;

#include <iostream>
#include <stdio.h>
#include <vector>


SimpleCloneTree::SimpleCloneTree (int numClones, std::vector<float> freqs) : m_numClones(numClones), m_vecNodes(2*numClones-1), m_vecLeafs(numClones) {
  // initialize tips
  for (int i=0; i<numClones; i++) {
    Clone *c = new Clone();
    c->label = i+1;
    c->freq = freqs[i];
    c->is_visible = true;
    m_vecNodes[i] = c;
    m_vecLeafs[i] = c;
  }
#ifdef DEBUG
  printNodes();
#endif
}

Clone* SimpleCloneTree::getRoot() {
  return m_root;
}

std::vector<Clone *> SimpleCloneTree::getLeafs() {
  return m_vecLeafs;
}

void SimpleCloneTree::generateRandomTopology(boost::function<float()>& random) {
  // generate N-1 internal nodes (each representing a coalescence event)
  int numClones = m_numClones;
  int k = numClones-1;
  int nextIndex = numClones;
  for (int i=0; i<numClones-1; i++) {
    // pick first random node (without replacement)
    int index1 = 2*i + random()*k--;
    Clone *p = m_vecNodes[index1];
    m_vecNodes[index1] = m_vecNodes[2*i];
    m_vecNodes[2*i] = p;
#ifdef DEBUG
    fprintf(stderr, "---\nIteration %d:\n", i);
    fprintf(stderr, "\tindex1: %d\n", index1);
    fprintf(stderr, "\tnode1: %d\n", p->label);
    printNodes();
#endif
    // pick second random node (without replacement)
    int index2 = 2*i+1 + random()*k;
    Clone *q = m_vecNodes[index2];
    m_vecNodes[index2] = m_vecNodes[2*i+1];
    m_vecNodes[2*i+1] = q;
#ifdef DEBUG
    fprintf(stderr, "\tindex2: %d\n", index2);
    fprintf(stderr, "\tnode2: %d\n", q->label);
    printNodes();
#endif
    // create new internal node
    Clone *n = new Clone();
    n->label = nextIndex+1;
    n->m_vecChildren.push_back(p);
    n->m_vecChildren.push_back(q);
    p->parent = n;
    q->parent = n;
    m_vecNodes[nextIndex++] = n;
#ifdef DEBUG
    fprintf(stderr, "\tnew internal node: %d\n", n->label);
    printNodes();
#endif
  }

  // generate a "healthy" clone as root node
  Clone *r = new Clone();
  r->label = 0;
  r->is_healthy = true;
  r->m_vecChildren.push_back(m_vecNodes[nextIndex-1]);
  r->parent = 0;
  m_vecNodes[nextIndex-1]->parent = r;
  m_vecNodes.push_back(r);
  m_root = r;
#ifdef DEBUG
  fprintf(stderr, "\twe have been ROOTed: %d\n", r->label);
  printNodes();
#endif
}

/** Place mutations randomly on the tree.
 * A set of mandatory mutations need to exist between clones,
 * otherwise they cannot be distinguished.
 */
void SimpleCloneTree::evolve(int numMutations, int numTransformingMutations, boost::function<float()>& random) {
  int num_randomMutations = numMutations - (m_numClones-1) - numTransformingMutations;
  putMandatoryMutations(numTransformingMutations, random);
  putRandomMutations(num_randomMutations, m_numClones, random);
}

void SimpleCloneTree::putMandatoryMutations(int num_transmuts, boost::function<float()>& random) {
fprintf(stderr, "\ndropping %d mandatory mutations (%d transforming + %d to assure distinct clones)...\n", m_numClones+num_transmuts-1, num_transmuts, m_numClones-1);
  int m=0;
  // MRCA of clones receives transforming mutations
  Clone *mrca = m_root->m_vecChildren[0];
  while (m<num_transmuts) {
fprintf(stderr, "\tdropping mutation %d on node %d\n", m, mrca->label);
    mrca->m_vecMutations.push_back(m);
    m++;
  }
  // each intermediate node needs to introduce a mutation
  _putMandatoryMutation(mrca, m, random);
fprintf(stderr, "%d mutations dropped.\n", m);
}

void SimpleCloneTree::_putMandatoryMutation(Clone *c, int& mutationId, boost::function<float()>& random) {
  Clone node = *c;
  if (node.isLeaf()) {
    return;
  }
  else {
    bool is_leftMutated = (random()<0.5); // throw a coin
    if (is_leftMutated) {
fprintf(stderr, "\tdropping mutation %d on node %d\n", mutationId, node.m_vecChildren[0]->label);
      node.m_vecChildren[0]->m_vecMutations.push_back(mutationId);
    }
    else {
fprintf(stderr, "\tdropping mutation %d on node %d\n", mutationId, node.m_vecChildren[1]->label);
      node.m_vecChildren[1]->m_vecMutations.push_back(mutationId);
    }
    mutationId++;
    _putMandatoryMutation(node.m_vecChildren[0], mutationId, random);
    _putMandatoryMutation(node.m_vecChildren[1], mutationId, random);
  }
}

void SimpleCloneTree::putRandomMutations(int num_mutations, int nextMutationId, boost::function<float()>& random) {
fprintf(stderr, "\nnow dropping %d random mutations...\n", num_mutations);
  while (num_mutations > 0) {
    // mutate random node
    int i = random()*(m_numClones+1);
    m_vecNodes[i]->m_vecMutations.push_back(nextMutationId);
fprintf(stderr, "\tdropping mutation %d on node %d\n", nextMutationId, m_vecNodes[i]->label);
    nextMutationId++;
    num_mutations--;
  }
}

void SimpleCloneTree::collapseZeroBranches(Clone *node) {
fprintf(stderr, "collapseZeroBranches(<Node(label=%d)>)\n", node->label);
  // check among children whether any has zero distance
  unsigned i=0;
  while ((i<node->m_vecChildren.size()) && ((int)(node->m_vecChildren[i]->distanceToParent())>0)) {
    i++;
  }
  if (i < node->m_vecChildren.size()) { // there is a child with zero distance
    Clone *child = node->m_vecChildren[i];
    // goodbye cruel world!
    child->replace(node);
    delete node;

fprintf(stderr, "<Node(label=%d)> (so fresh & so clean)\n", child->label);
fprintf(stderr, "\tparent <Node(label=%d)>\n", child->parent->label);
for (unsigned k=0; k<child->m_vecChildren.size(); k++) {
fprintf(stderr, "\tchild <Node(label=%d)>\n", child->m_vecChildren[k]->label);
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

// use me for debugging :-)
void SimpleCloneTree::printNodes() {
  for (unsigned i=0; i<m_vecNodes.size(); i++) { fprintf(stderr, "|%2u ", i); }; fprintf(stderr, "|\n");
  for (unsigned i=0; i<m_vecNodes.size(); i++) { fprintf(stderr, "+---"); }; fprintf(stderr, "+\n");
  for (unsigned i=0; i<m_vecNodes.size(); i++) {
    if (m_vecNodes[i] > 0) { fprintf(stderr, "|%2d ", m_vecNodes[i]->label); }
    else { fprintf(stderr, "| - "); }
  };
  fprintf(stderr, "|\n");
}

void SimpleCloneTree::printDot(Clone *node, std::ostream& os) {
  os << "digraph G {\n";
  _printDotRecursive(node, os);
  os << "}\n";
}

void SimpleCloneTree::_printDotRecursive(Clone *node, std::ostream& os) {
  if (node->is_healthy) {
    os << "\t" << node->label << " [style=filled,color=limegreen];" << std::endl;
  } else if (node->is_visible) {
    os << "\t" << node->label << " [style=filled,color=tomato];" << std::endl;
  }
  for (unsigned i=0; i<node->m_vecChildren.size(); ++i) {
    Clone *child = node->m_vecChildren[i];
    float edgeWeight = child->distanceToParent();
    os << "\t" << node->label << " -> " << child->label;
    if (edgeWeight > 0.0) {
      os << "[style=bold,label=" << edgeWeight << "]";
    }
    os << ";" << std::endl;

    _printDotRecursive(child, os);
  }
}

/** Generates the string representation of a (sub)tree in Newick notation. */
void SimpleCloneTree::printNewick(Clone *node, std::ostream& os) {
  _printNewickRecursive(node, true, os);
  os << ";";
}

/** Prints the string representation of the tree in Newick notation. (internal) */
void SimpleCloneTree::_printNewickRecursive(Clone *node, bool isFirstChild, std::ostream& os) {
  if (!isFirstChild) { os << ","; }
  if (node->getChildren().size() > 0) {
    os << "(";
    isFirstChild = true;
    for (unsigned i=0; i<node->getChildren().size(); i++) {
      _printNewickRecursive(node->getChildren()[i], isFirstChild, os);
      isFirstChild = false;
    }
    os << ")" << node->label << ":" << node->distanceToParent()+0.05;
  }
  else {
    os << node->label << ":" << node->distanceToParent()+0.05;
  }
}
