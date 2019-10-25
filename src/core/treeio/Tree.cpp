#include "Tree.hpp"
#include "../stringio.hpp"
using stringio::format;
#include <boost/lexical_cast.hpp>
#include <fstream>

using namespace std;

namespace treeio {

template<typename TNodeType>
Tree<TNodeType>::Tree() : m_numNodes(0), m_numVisibleNodes(0), m_vecNodes(0) {}

template<typename TNodeType>
Tree<TNodeType>::Tree(int n_nodes) :
  m_numNodes(n_nodes),
  m_numVisibleNodes(n_nodes),
  m_vecNodes(n_nodes),
  m_numMutations(0)
{
  // initialize nodes
  for (int i=0; i<n_nodes; i++) {
    shared_ptr<TNodeType> n(new TNodeType());
    n->index = i;
    n->label = boost::lexical_cast<std::string>(i+1);
    n->is_visible = true;
    n->length = 1;
    n->weight = 0.0;
    m_vecNodes[i] = n;
  }
#ifdef DEBUG
  _printNodes();
#endif
}

template<typename TNodeType>
Tree<TNodeType>::Tree(const parse::node& root) :
  m_numNodes(0),
  m_numVisibleNodes(0),
  m_numMutations(0)
{
  std::shared_ptr<TNodeType> root_node = _adaptNode(root);
  root_node->parent = 0;
  this->m_root = root_node;
}

template<typename TNodeType>
Tree<TNodeType>::Tree(const string filename) :
  m_numNodes(0),
  m_numVisibleNodes(0),
  m_numMutations(0)
{
  fprintf(stderr, "\nReading tree from file '%s'...\n", filename.c_str());
  parse::node root;
  readNewick(filename, root);
  std::shared_ptr<TNodeType> root_node = _adaptNode(root);
  root_node->parent = 0;

  this->m_root = root_node;
  _printTreeInfo();
}

template<typename T>
Tree<T>::~Tree() {
  // NOTE:this loop causes a segfault in test bamio/bulk
  /*for (unsigned i=0; i<m_vecNodes.size(); ++i) {
    delete(m_vecNodes[i]);
  }*/
}

template<typename TNodeType>
shared_ptr<TNodeType> Tree<TNodeType>::_adaptNode(const parse::node& node) {
  //T *n = new T();
  shared_ptr<TNodeType> n(new TNodeType());
  n->index = m_numNodes++;
  n->label = node.label;
  n->length = node.length;
  n->is_visible = (node.label.size() > 0); // unlabeled nodes are invisible
  m_numVisibleNodes += n->is_visible ? 1 : 0;
  this->m_vecNodes.push_back(n);

  for (unsigned i=0; i<node.children.size(); ++i) {
    shared_ptr<TNodeType> child_node = _adaptNode(node.children[i]);
    n->m_vecChildren.push_back(child_node);
    child_node->parent = n;
  }

  return n;
}

template<typename TNodeType>
double Tree<TNodeType>::getTotalBranchLength() {
  double branch_len = 0;
  for (unsigned i=0; i<m_vecNodes.size(); ++i) {
    branch_len += m_vecNodes[i]->length;
  }
  return branch_len;
}

template<typename TNodeType>
vector<double> Tree<TNodeType>::getAbsoluteBranchLengths() {
  vector<double> vec_len = vector<double>(m_vecNodes.size());
  for (unsigned i=0; i<this->m_vecNodes.size(); ++i) {
    vec_len[i] = this->m_vecNodes[i]->length;
  }
  return vec_len;
}

template<typename TNodeType>
vector<double> Tree<TNodeType>::getRelativeBranchLengths() {
  double total_branch_len = this->getTotalBranchLength();
  vector<double> branch_lens = vector<double>(m_vecNodes.size());
  for (unsigned i=0; i<this->m_vecNodes.size(); ++i) {
    branch_lens[i] = this->m_vecNodes[i]->length / total_branch_len;
  }
  return branch_lens;
}

template<typename TNodeType>
vector<shared_ptr<TNodeType>> Tree<TNodeType>::getVisibleNodes() {
  vector<shared_ptr<TNodeType>> vis_nodes;
  for (unsigned i=0; i<m_vecNodes.size(); ++i) {
    shared_ptr<TNodeType> node = m_vecNodes[i];
    //if (node->is_visible && (node!=m_root)) { // root represents contaminant healthy cells
    if (node->is_visible) {
      vis_nodes.push_back(node);
    }
  }
  return vis_nodes;
}

template<typename TNodeType>
vector<int> Tree<TNodeType>::getVisibleNodesIdx() {
  vector<int> vis_nodes_idx;
  for (auto node : m_vecNodes) {
    if (node->is_visible)
      vis_nodes_idx.push_back(node->index);
  }
  return vis_nodes_idx;
}

template<typename TNodeType>
vector<shared_ptr<TNodeType>> Tree<TNodeType>::getNodesPreOrder() {
  vector<shared_ptr<TNodeType>> nodes;
  getNodesPreOrderRec(this->m_root, nodes);
  return nodes;
}

template<typename TNodeType>
void Tree<TNodeType>::getNodesPreOrderRec (
  const shared_ptr<TNodeType> n,
  vector<shared_ptr<TNodeType>>& nodes)
{
  nodes.push_back(n);
  for (auto const c : n->m_vecChildren) {
    getNodesPreOrderRec(c, nodes);
  }
}


template<typename TNodeType>
void Tree<TNodeType>::generateRandomTopology(function<double()>& rng) {
  // TODO: call appropriate method based on user params
  if (true) {
    generateRandomTopologyInternalNodes(rng);
  }
  else {
    generateRandomTopologyLeafsOnly(rng);
  }
}

template<typename TNodeType>
void Tree<TNodeType>::generateRandomTopologyInternalNodes(function<double()>& random) {
  // create root node
  shared_ptr<TNodeType> r(new TNodeType());
  r->label = "0";
  //r->is_healthy = true;
  r->parent = 0;
  r->index = m_vecNodes.size();
  r->is_visible = true;
  m_vecNodes.push_back(r);
  m_numNodes++;
  m_numVisibleNodes++;
  m_root = r;
  // first clone becomes child of root
  Clone::setParent(m_vecNodes[0], r);

  // pick parent for each clone
  vector<shared_ptr<TNodeType>> parents;
  parents.push_back(m_vecNodes[0]);
  for (int i=1; i<m_numNodes-1; ++i) {
    shared_ptr<TNodeType> n = m_vecNodes[i];
    int p_index = random()*i;
    shared_ptr<TNodeType> p = parents[p_index];
// cerr << *n << " gets parent " << *p << std::endl;
    Clone::setParent(n, p);
    parents.push_back(n);
  }
}

template<typename TNodeType>
void Tree<TNodeType>::generateRandomTopologyLeafsOnly(function<double()>& random) {
  // generate N-1 internal nodes (each representing a coalescence event)
  int numNodes = m_numVisibleNodes;
  int k = numNodes-1;
  int nextIndex = numNodes;
  for (int i=0; i<numNodes-1; i++) {
    // pick first random node (without replacement)
    int index1 = 2*i + random()*k--;
    shared_ptr<TNodeType> p = m_vecNodes[index1];
    m_vecNodes[index1] = m_vecNodes[2*i];
    m_vecNodes[2*i] = p;
    p->index = 2*i;
#ifdef DEBUG
    fprintf(stderr, "---\nIteration %d:\n", i);
    fprintf(stderr, "\tindex1: %d\n", index1);
    fprintf(stderr, "\tnode1: %s\n", p->label.c_str());
    _printNodes();
#endif
    // pick second random node (without replacement)
    int index2 = 2*i+1 + random()*k;
    shared_ptr<TNodeType> q = m_vecNodes[index2];
    m_vecNodes[index2] = m_vecNodes[2*i+1];
    m_vecNodes[2*i+1] = q;
    q->index = 2*i+1;
#ifdef DEBUG
    fprintf(stderr, "\tindex2: %d\n", index2);
    fprintf(stderr, "\tnode2: %s\n", q->label.c_str());
    _printNodes();
#endif
    // create new internal node
    shared_ptr<TNodeType> n(new TNodeType());
    n->index = nextIndex++;
    n->label = boost::lexical_cast<string>(n->index+1);
    n->length = 1;
    n->weight = 0.0;
    n->m_vecChildren.push_back(p);
    n->m_vecChildren.push_back(q);
    p->parent = n;
    q->parent = n;
    m_vecNodes.push_back(n);
    m_numNodes++;
#ifdef DEBUG
    fprintf(stderr, "\tnew internal node: %s\n", n->label.c_str());
    _printNodes();
#endif
  }

  // generate a "healthy" clone as root node
  shared_ptr<TNodeType> r(new TNodeType());
  r->index = nextIndex;
  r->label = "0";
  //r->is_healthy = true;
  r->is_visible = true;
  r->m_vecChildren.push_back(m_vecNodes[nextIndex-1]);
  r->parent = 0;
  m_vecNodes[nextIndex-1]->parent = r;
  m_vecNodes.push_back(r);
  m_numNodes++;
  m_numVisibleNodes++;
  m_root = r;
#ifdef DEBUG
  fprintf(stderr, "\twe have been ROOTed: %s\n", r->label.c_str());
  _printNodes();
#endif
}

/** Shrink/expand branch length by a random factor */
template <typename TNodeType>
void Tree<TNodeType>::varyBranchLengths(
  function<double()>& random_double
)
{
  shared_ptr<TNodeType> root = this->m_root;
  _varyBranchLengthsRec(root, random_double);
  // normalize branch lengths to MRCA's length
  /*double len_mrca = this->m_root->m_vecChildren[0]->length;
  for (T* node : this->m_vecNodes) {
    node->length /= len_mrca;
  }*/
}

template <typename TNodeType>
void Tree<TNodeType>::_varyBranchLengthsRec(
  shared_ptr<TNodeType> node,
  function<double()>& random_double
)
{
  for (auto child : node->m_vecChildren) {
    double scale_factor = random_double();
    child->length *= scale_factor;
    _varyBranchLengthsRec(child, random_double);
  }
  // normalize branch lengths to MRCA's length
  //double len_mrca = this->m_root->m_vecChildren[0]->length;
  //for (T* node : this->m_vecNodes) {
  //  node->length /= len_mrca;
  //}
}

/** Assign random weights to visible nodes */
template<typename TNodeType>
void Tree<TNodeType>::assignWeights(vector<double> w) {
  string str_w = format("%.4f", w[0]);
  for_each(w.begin()+1, w.end(), [&] (double p) {
    str_w += format(", %.4f", p);
  });
  fprintf(stderr, "Assigning weights:\n\t[ %s ]\n", str_w.c_str());

  // note: first weight assigned to root node -> normal cell contamination
  this->m_root->weight = w[0];
  long num_weights = w.size();
  if (this->m_numVisibleNodes != num_weights) {
    fprintf(stderr, "[ERROR] number of node weights (%ld) must equal number of visible nodes (%d)\n", num_weights, m_numVisibleNodes);
  }
  auto it_w = w.begin()+1;
  for (auto child : this->m_root->m_vecChildren) {
    _assignWeightsRec(child, it_w, w);
  }
/*  long i = 1;
  for (auto node : m_vecNodes) {
    if (node->is_visible) {
      if (i == num_weights) {
        fprintf(stderr, "[ERROR] number of visible nodes exceeds number of weights (%ld)\n", num_weights);
      }
      node->weight = w[i++];
    }
  }
  if (i < num_weights) {
    fprintf(stderr, "[WARN] number of visible nodes (%ld) less than number of weights (%ld)\n", i, num_weights);
  }*/
}

/** Assign weights recursively (pre-order traversal) */
template<typename TNodeType>
void Tree<TNodeType>::_assignWeightsRec(
  shared_ptr<TNodeType> node,
  vector<double>::iterator &it_w,
  vector<double> w)
{
  if (node->is_visible) {
    if (it_w == w.end()) {
      fprintf(stderr, "[ERROR] end of node weights reached before end of nodes.\n");
      return;
    }
    node->weight = *it_w;
    it_w++;
  }

  for (auto child : node->m_vecChildren) {
    _assignWeightsRec(child, it_w, w);
  }
}

/** Place mutations randomly on the tree.
 * Relative branch lengths are used as prior probabilities during assignment,
 * so longer branches receive more mutations, compared to shorter ones.
 *
 * After dropping mutations, branch lengths are reassigned as number of mutations.
 */
template<typename TNodeType>
void Tree<TNodeType>::dropSomaticMutations(
  int n_mutations,
  int n_transforming,
  RandomNumberGenerator &rng)
{
  this->m_numMutations = n_mutations;
#ifdef DEBUG
  fprintf(stderr, "Dropping %d mutations (%d transforming)...\n", n_mutations, n_transforming);
#endif
  dropTransformingMutations(n_transforming);
  int next_mut_id = n_transforming; // identifier for next mutation
#ifdef DEBUG
  //fprintf(stderr, "Now dropping mandatory mutations...\n");
#endif
  /* TODO: this makes it difficult to have internal clones
    (represented by zero-branch from "dummy" internal node) */
  // each clone needs at least one private mutation
  //dropMandatoryMutations(m_root, next_mut_id);
  // remaining mutations are dropped randomly, proportional to branch_length
  //int n_random = n_mutations - n_transforming - m_numVisibleNodes;
  int n_random = n_mutations - n_transforming;
#ifdef DEBUG
  fprintf(stderr, "Now dropping %d random mutations...\n", n_random);
#endif
  dropRandomMutations(n_random, ++next_mut_id, rng);

  // relabel mutations to pre-order sequence
  next_mut_id = 0;
  _relabelMutationsRec(m_root, next_mut_id);

  // reset branch lengths to #mutations
  for (auto node : m_vecNodes) {
    node->length = node->m_vec_mutations.size();
  }
}

/** Drop transforming mutations on root (i.e. healthy) node. 
 */
template<typename TNodeType>
void Tree<TNodeType>::dropTransformingMutations(int n_mutations) {
  shared_ptr<TNodeType> node = m_root;
//  vector<shared_ptr<TNodeType>> topNodes = m_root->m_vecChildren;
//  assert ( topNodes.size() == 1 && "Root should have only one child." );
//  for (unsigned i=0; i<topNodes.size(); ++i) {
//    shared_ptr<TNodeType> node = topNodes[i];
    for (int m=0; m<n_mutations; ++m) {
cerr << "\tDropping mutation " << m << " on " << *node << endl;;
//fprintf(stderr, "\tDropping mutation %d on Clone<label=%s>\n", m, c->label.c_str());
      node->m_vec_mutations.push_back(m);
    }
//  }
}

/** Drop one mutation on a given node and repeat for children. */
template<typename TNodeType>
void Tree<TNodeType>::dropMandatoryMutations(shared_ptr<TNodeType> node, int &mutation_id) {
  if (node!=m_root) {
cerr << "\tDropping mutation " << mutation_id << " on " << *node << endl;
//fprintf(stderr, "\tDropping mutation %d on Clone<label=%d>\n", mutationId, clone->label);
    node->m_vec_mutations.push_back(mutation_id);
  }
  else {
    mutation_id--;
  }
  for (unsigned i=0; i<node->m_vecChildren.size(); ++i) {
    dropMandatoryMutations(node->m_vecChildren[i], ++mutation_id);
  }
}

/** Drop mutations randomly along tree. */
template<typename TNodeType>
void 
Tree<TNodeType>::dropRandomMutations (
  int n_mutations, 
  int &mutation_id, 
  RandomNumberGenerator& rng
) {
  long n_nodes = this->m_vecNodes.size();
  function<int()> f_random_index;
  // are branch lengths specified?
  double tot_len = this->getTotalBranchLength();
  vector<double> vec_branch_len = this->getAbsoluteBranchLengths();
  if (tot_len == 0) { // if tree has no branch length, assign equal weight to each
    vec_branch_len = vector<double>(n_nodes, 1);
  }
  // avoid MRCA node receiving random mutations
  shared_ptr<TNodeType> mrca = this->m_root;
  vec_branch_len[mrca->index] = 0.0;
  f_random_index = rng.getRandomIndexWeighted(vec_branch_len);
  // drop mutations randomly, but in proportion to branch length
  for (int i=0; i<n_mutations; ++i) {
    // pick clone to mutate
    long idx_node = f_random_index();
    shared_ptr<TNodeType> node = m_vecNodes[idx_node];
//cerr << "\tDropping mutation " << mutation_id << " on " << *node << endl;
//fprintf(stderr, "\tDropping mutation %d on Clone<label=%d>\n", mutationId, c->label);
    node->m_vec_mutations.push_back(mutation_id++);
  }
}

/** Reset mutation ids to follow pre-oder traversal sequence. */
template<typename TNodeType>
void Tree<TNodeType>::_relabelMutationsRec(shared_ptr<TNodeType> node, int &next_mut_id) {
  // relabel this node's mutations
  for (int &m_id : node->m_vec_mutations)
    m_id = next_mut_id++;
  // relabel children's mutations
  for (shared_ptr<TNodeType> child : node->m_vecChildren)
    _relabelMutationsRec(child, next_mut_id);
}

template<typename TNodeType>
void Tree<TNodeType>::printNewick(const string filename) {
  ofstream fs;
  fs.open(filename);
  printNewick(fs);
  fs.close();
}

template<typename TNodeType>
void Tree<TNodeType>::printNewick(ostream& os) {
  printNewick(this->m_root, os);
  os << ';' << endl;
}

template<typename TNodeType>
void Tree<TNodeType>::printNewick(shared_ptr<TNodeType> node, ostream& os, bool first) {
  bool has_children = (node->m_vecChildren.size() > 0);
  if (!first) {
    os << ","; }
  if (has_children) {
    os << "(";
    bool first_child = true;
    for (auto child : node->m_vecChildren) {
      printNewick(child, os, first_child);
      first_child = false;
    }
    os << ")";
  }
  os << format("%s:%g", node->label.c_str(), node->length);
}

template<typename TNodeType>
void Tree<TNodeType>::printNexus(const string filename) {
  ofstream fs;
  fs.open(filename);
  printNexus(fs);
  fs.close();
}

template<typename TNodeType>
void Tree<TNodeType>::printNexus(ostream& os) {
  os << "#NEXUS" << endl;
  os << "BEGIN DATA;" << endl;
  os << "  Dimensions" 
     << " ntax=" << m_numNodes 
     << " nchar=" << m_numMutations 
     << ';' << endl;
  os << "  Format datatype=standard;" << endl;
  printNexusMatrix( os );
  os << "END;" << endl;
  os << "BEGIN TREES;" << endl;
  // printNexusTranslate( os );
  os << "  Tree clones = ";
  printNexusTree( this->m_root, os );
  os << ";" << endl;
  os << "END;" << endl;
}

template<typename TNodeType>
void Tree<TNodeType>::printNexusTranslate (
  ostream& os
)
{
  auto node = this->m_vecNodes[0];
  os << "  Translate" << endl << (
    node->label.length() > 0 ?
    format("    %d %s", node->index, node->label.c_str()) :
    format("    %d internal_%d", node->index, node->index)
    );
  for ( size_t i=1; i<this->m_vecNodes.size(); i++ ) {
    node = this->m_vecNodes[i];
    os << ',' << endl << (
      node->label.length() > 0 ?
      format("    %d %s", node->index, node->label.c_str()) :
      format("    %d internal_%d", node->index, node->index)
    );
  }
  os << endl << "  ;" << endl;
}

template<typename TNodeType>
void Tree<TNodeType>::printNexusMatrix (
  ostream& os
)
{
  os << "  Matrix" << endl;
  os << "    ['1's indicate on which branch a mutation *first* occurred]" << endl;
  for ( size_t i=0; i<m_vecNodes.size(); i++ ) {
    auto node = m_vecNodes[i];
    string lbl = (
      node->label.length() > 0 ?
      format("    %s", node->label.c_str()) :
      format("    I%d", node->index)
    );
    os << "    " << lbl;
    // add padding to print matrix horizontally aligned
    for ( size_t j=0; j<(10-lbl.length()); j++) {
      os << ' ';
    }
    // set matrix line with new mutations for this node
    vector<int> mm(m_numMutations, 0);
    for ( auto const k : node->m_vec_mutations ) {
      mm[k] = 1;
    }
    // print mutation matrix as 0/1 values
    for ( size_t l=0; l<mm.size(); l++) {
      os << mm[l];
    }
    os << endl;
  }
  os << endl << "  ;" << endl;
}

template<typename TNodeType>
void Tree<TNodeType>::printNexusTree (
  shared_ptr<TNodeType> node, 
  ostream& os, 
  bool first
)
{
  bool has_children = (node->m_vecChildren.size() > 0);
  if (!first) {
    os << ","; }
  if (has_children) {
    os << "(";
    bool first_child = true;
    for (auto child : node->m_vecChildren) {
      printNexusTree(child, os, first_child);
      first_child = false;
    }
    os << ")";
  }
  string lbl = (
    node->label.length() > 0 ? 
    node->label : 
    format("I%d", node->index)
  );
  os << format("%s:%g", lbl.c_str(), node->length);
}

/** Print tree a graph in DOT format. */
template<typename TNodeType>
void Tree<TNodeType>::printDot(const string filename) {
  ofstream fs;
  fs.open(filename);
  printDot(this->m_root, fs);
  fs.close();
}

/** Print tree a graph in DOT format. */
template<typename TNodeType>
void Tree<TNodeType>::printDot(shared_ptr<TNodeType> node, std::ostream& os) {
  os << "digraph G {\n";
  _printDotRec(node, os);
  os << "}\n";
}

template<typename TNodeType>
void Tree<TNodeType>::_printDotRec(shared_ptr<TNodeType> node, std::ostream& os) {
  //auto node_lbl = boost::format("\"%s\\ni:%d\\nw:%.4f\"") % node->label % node->index % node->weight;
  string node_lbl = format("\"%s\\ni:%d\"", node->label.c_str(), node->index);
  os << "\t" << node->index << "[label=" << node_lbl;
  if (node==this->m_root) {
    os << ",style=filled,color=limegreen";
  } else if (node->is_visible) {
    os << ",style=filled,color=tomato";
  }
  os << "];" << std::endl;
  for (auto child : node->m_vecChildren) {
    float edgeMut = child->distanceToParent();
    float edgeLen = child->length;
    os << "\t" << node->index << " -> " << child->index;
    if (edgeLen > 0.0) {
      string lbl = format("%0.2f (%0.0f)", edgeLen, edgeMut);
      os << "[style=bold,label=<<font point-size=\"10\">" << lbl << "</font>>]";
    }
    os << ";" << std::endl;

    _printDotRec(child, os);
  }
}

/** outputs boolean matrix of mutational states for visible clones */
template<>
void Tree<Clone>::writeMutationMatrix(ostream& os) {
  // setup matrix
  int num_nodes = this->m_numNodes;
  int num_mutations = this->m_numMutations;
  vector<vector<bool>> mm(num_nodes, vector<bool>(num_mutations, false));
  this->m_root->populateMutationMatrixRec(mm);
  // write matrix
  vector<shared_ptr<Clone>> vec_vis_clones = this->getVisibleNodes();
  for (auto clone : vec_vis_clones) {
    os << clone->label;
    for (auto cell : mm[clone->index])
      os << "," << cell;
    os << endl;
  }
}

/** outputs boolean matrix of mutational states for visible clones */
template<>
void Tree<Clone>::writeMutationMatrix(const string filename) {
  ofstream fs;
  fs.open(filename);
  writeMutationMatrix(fs);
  fs.close();
}

// use me for debugging :-)
template<typename T>
void Tree<T>::_printNodes() {
  for (unsigned i=0; i<m_vecNodes.size(); i++) { fprintf(stderr, "|%2u ", i); }; fprintf(stderr, "|\n");
  for (unsigned i=0; i<m_vecNodes.size(); i++) { fprintf(stderr, "+---"); }; fprintf(stderr, "+\n");
  for (unsigned i=0; i<m_vecNodes.size(); i++) {
    if (m_vecNodes[i]) { fprintf(stderr, "|%s ", m_vecNodes[i]->label.c_str()); }
    else { fprintf(stderr, "| - "); }
  };
  fprintf(stderr, "|\n");
}

template<typename TNodeType>
void Tree<TNodeType>::_printTreeInfo() {
  fprintf(stderr, "\n");
  fprintf(stderr, "m_numNodes: %d\n", this->m_numNodes);
  fprintf(stderr, "m_numVisibleNodes: %d\n", this->m_numVisibleNodes);
  fprintf(stderr, "\n");
}

} // namespace treeio

// instantiate usable classes so they can be picked up by the linker
template class treeio::Tree<Clone>;