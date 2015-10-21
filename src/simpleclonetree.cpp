#include "simpleclonetree.hpp"
//#include "node.hpp"
#include <boost/random.hpp>
// Choosing the random number generator. (mt19937: Mersenne-Twister)
typedef boost::mt19937 base_generator_type;

#include <iostream>
#include <stdio.h>
#include <vector>

SimpleCloneTree::SimpleCloneTree (int num_clones, std::vector<float> freqs) : Tree(num_clones) {
  // initialize tips
  for (int i=0; i<num_clones; i++) {
    Clone *c = new Clone();
    c->label = i+1;
    c->freq = freqs[i];
    c->is_healthy = false;
    m_vecNodes[i] = c;
  }

  printNodes();
}

void SimpleCloneTree::generateRandomTopology(long seed) {
  fprintf(stderr, "seed: %ld\n", seed);
  // init random number generator
  boost::random::mt19937 rng(seed);
  boost::uniform_real<> uni_dist(0, 1);
  boost::variate_generator<base_generator_type&, boost::uniform_real<> > uni(rng, uni_dist);
  // take that baby for a spin
  //for (int i=0; i<10; i++) { fprintf(stderr, "%d\n", m_vecNodes[uni()*m_numClones].label); }

  // generate N-1 internal nodes (each representing a coalescence event)
  int numClones = m_numTips;
  int k = numClones-1;
  int nextIndex = numClones;
  for (int i=0; i<numClones-1; i++) {
fprintf(stderr, "---\nIteration %d:\n", i);
    // pick first random node (without replacement)
    int index1 = 2*i + uni()*k--;
fprintf(stderr, "\tindex1: %d\n", index1);
    Node *p = m_vecNodes[index1];
fprintf(stderr, "\tnode1: %d\n", p->label);
    m_vecNodes[index1] = m_vecNodes[2*i];
    m_vecNodes[2*i] = p;
this->printNodes();
    // pick second random node (without replacement)
    int index2 = 2*i+1 + uni()*k;
fprintf(stderr, "\tindex2: %d\n", index2);
    Node *q = m_vecNodes[index2];
fprintf(stderr, "\tnode2: %d\n", q->label);
    m_vecNodes[index2] = m_vecNodes[2*i+1];
    m_vecNodes[2*i+1] = q;
this->printNodes();
    // create new internal node
    Node *n = new Node();
    n->label = nextIndex+1;
    n->m_vecChildren.push_back(p);
    n->m_vecChildren.push_back(q);
    m_vecNodes[nextIndex++] = n;
fprintf(stderr, "\tnew internal node: %d\n", n->label);
this->printNodes();
  }

  // generate a "healthy" clone as root node
  Clone *r = new Clone();
  r->label = 0;
  r->is_healthy = true;
  r->m_vecChildren.push_back(m_vecNodes[nextIndex-1]);
  m_vecNodes.push_back(r);
  m_root = r;
fprintf(stderr, "\twe have been ROOTed: %d\n", r->label);
this->printNodes();
}
