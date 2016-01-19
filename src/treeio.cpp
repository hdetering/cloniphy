#define DEBUG
#include "clone.hpp"
#include "treeio.hpp"
#include <boost/lexical_cast.hpp>
#include <boost/spirit/include/qi.hpp>
#include <boost/fusion/include/adapt_struct.hpp>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <streambuf>
#include <string>
#include <vector>

using namespace std;

Node::Node() {
  parent = NULL;
}

Node::~Node() {}

// streaming operator for easy printing
ostream& operator<<(std::ostream& stream, const Node& node) {
  stream << "<Node(index=" << node.index << ", label='" << node.label << "', length=" << node.length << ")>";
  return stream;
}

bool Node::isRoot() {
  return (parent == NULL);
}

template<typename T>
Tree<T>::Tree() : m_numNodes(0), m_vecNodes(0) {}

template<typename T>
Tree<T>::Tree(int n_nodes) : m_numNodes(n_nodes), m_vecNodes(n_nodes) {
  // initialize nodes
  for (int i=0; i<n_nodes; i++) {
    T *n = new T();
    n->index = i+1;
    n->label = boost::lexical_cast<std::string>(i+1);
    n->is_visible = true;
    m_vecNodes[i] = n;
  }
#ifdef DEBUG
  _printNodes();
#endif
}

template<typename T>
Tree<T>::Tree(const treeio::node& root) {
  this->m_numNodes = 0;
  T *root_node = _adaptNode(root);
  root_node->parent = 0;
  this->m_root = root_node;
}

template<typename T>
Tree<T>::~Tree() {
  for (unsigned i=0; i<m_vecNodes.size(); ++i) {
    delete(m_vecNodes[i]);
  }
}

template<typename T>
T* Tree<T>::_adaptNode(const treeio::node& node) {
  T *n = new T();
  n->index = m_numNodes++;
  n->label = node.label;
  n->length = node.length;
  n->is_visible = (node.label.size() > 0); // unlabeled nodes are invisible
  this->m_vecNodes.push_back(n);

  for (unsigned i=0; i<node.children.size(); ++i) {
    T *child_node = _adaptNode(node.children[i]);
    n->m_vecChildren.push_back(child_node);
    child_node->parent = n;
  }

  return n;
}

template<typename T>
vector<T *> Tree<T>::getVisibleNodes() {
  vector<T *> vis_nodes;
  for (unsigned i=0; i<m_vecNodes.size(); ++i) {
    T *node = m_vecNodes[i];
    if (node->is_visible && (node!=m_root)) {
      vis_nodes.push_back(node);
    }
  }
  return vis_nodes;
}

template<typename T>
void Tree<T>::generateRandomTopology(boost::function<float()>& rng) {
  // TODO: call appropriate method based on user params
  if (true) {
    generateRandomTopologyInternalNodes(rng);
  }
  else {
    generateRandomTopologyLeafsOnly(rng);
  }
}

template<typename T>
void Tree<T>::generateRandomTopologyInternalNodes(boost::function<float()>& random) {
  // create root node
  T *r = new T();
  r->label = "0";
  //r->is_healthy = true;
  r->parent = NULL;
  m_vecNodes.push_back(r);
  m_root = r;
  // first clone becomes child of root
  m_vecNodes[0]->setParent(r);

  // pick parent for each clone
  std::vector<T*> parents;
  parents.push_back(m_vecNodes[0]);
  for (int i=1; i<m_numNodes; ++i) {
    T *n = m_vecNodes[i];
    int p_index = random()*i;
    T *p = parents[p_index];
std::cerr << *n << " gets parent " << *p << std::endl;
    n->setParent(p);
    parents.push_back(n);
  }
}

template<typename T>
void Tree<T>::generateRandomTopologyLeafsOnly(boost::function<float()>& random) {
  // generate N-1 internal nodes (each representing a coalescence event)
  int numNodes = m_numNodes;
  int k = numNodes-1;
  int nextIndex = numNodes;
  for (int i=0; i<numNodes-1; i++) {
    // pick first random node (without replacement)
    int index1 = 2*i + random()*k--;
    T *p = m_vecNodes[index1];
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
    T *q = m_vecNodes[index2];
    m_vecNodes[index2] = m_vecNodes[2*i+1];
    m_vecNodes[2*i+1] = q;
#ifdef DEBUG
    fprintf(stderr, "\tindex2: %d\n", index2);
    fprintf(stderr, "\tnode2: %s\n", q->label.c_str());
    _printNodes();
#endif
    // create new internal node
    T *n = new T();
    n->label = boost::lexical_cast<string>(++nextIndex);
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
  T *r = new T();
  r->label = "0";
  //r->is_healthy = true;
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
template<typename T>
void Tree<T>::evolve(int n_mutations, int n_transforming, boost::function<float()>& random) {
#ifdef DEBUG
  fprintf(stderr, "Dropping %d mutations (%d transforming)...\n", n_mutations, n_transforming);
#endif
  dropTransformingMutations(n_transforming);
  int next_mut_id = n_transforming; // identifier for next mutation
#ifdef DEBUG
  fprintf(stderr, "Now dropping mandatory mutations...\n");
#endif
  // each clone needs at least one private mutation
  dropMandatoryMutations(m_root, next_mut_id);
  // remaining mutations are dropped randomly
  int n_random = n_mutations - n_transforming - m_numNodes;
#ifdef DEBUG
  fprintf(stderr, "Now dropping %d random mutations...\n", n_random);
#endif
  dropRandomMutations(n_random, ++next_mut_id, random);
}

/** Drop transforming mutations on immediate children of root node. */
template<typename T>
void Tree<T>::dropTransformingMutations(int n_mutations) {
  vector<T *> topNodes = m_root->m_vecChildren;
  for (unsigned i=0; i<topNodes.size(); ++i) {
    T *node = topNodes[i];
    for (int m=0; m<n_mutations; ++m) {
cerr << "\tDropping mutation " << m << " on " << *node << endl;;
//fprintf(stderr, "\tDropping mutation %d on Clone<label=%s>\n", m, c->label.c_str());
      node->m_vecMutations.push_back(m);
    }
  }
}

/** Drop one mutation on a given node and repeat for children. */
template<typename T>
void Tree<T>::dropMandatoryMutations(T *node, int &mutation_id) {
  if (node!=m_root) {
cerr << "\tDropping mutation " << mutation_id << " on " << *node << endl;
//fprintf(stderr, "\tDropping mutation %d on Clone<label=%d>\n", mutationId, clone->label);
    node->m_vecMutations.push_back(mutation_id);
  }
  else {
    mutation_id--;
  }
  for (unsigned i=0; i<node->m_vecChildren.size(); ++i) {
    dropMandatoryMutations(node->m_vecChildren[i], ++mutation_id);
  }
}

/** Drop mutations randomly along tree. */
template<typename T>
void Tree<T>::dropRandomMutations(int n_mutations, int &mutation_id, boost::function<float()>& random) {
  for (int i=0; i<n_mutations; ++i) {
    // pick random clone to mutate
    T *node = m_vecNodes[random()*m_numNodes];
cerr << "\tDropping mutation " << mutation_id << " on " << *node << endl;
//fprintf(stderr, "\tDropping mutation %d on Clone<label=%d>\n", mutationId, c->label);
    node->m_vecMutations.push_back(mutation_id++);
  }
}

/** Print tree a graph in DOT format. */
template<typename T>
void Tree<T>::printDot(T *node, std::ostream& os) {
  os << "digraph G {\n";
  _printDotRecursive(node, os);
  os << "}\n";
}

template<typename T>
void Tree<T>::_printDotRecursive(T *node, std::ostream& os) {
  if (node->isRoot()) {
    os << "\t" << node->index << " [style=filled,color=limegreen];" << std::endl;
  } else if (node->is_visible) {
    os << "\t" << node->index << " [style=filled,color=tomato];" << std::endl;
  }
  for (unsigned i=0; i<node->m_vecChildren.size(); ++i) {
    T *child = node->m_vecChildren[i];
    float edgeWeight = child->distanceToParent();
    os << "\t" << node->index << " -> " << child->index;
    if (edgeWeight > 0.0) {
      os << "[style=bold,label=" << edgeWeight << "]";
    }
    os << ";" << std::endl;
    os << "\t" << node->index << " [label=\"" << node->label << "\"];" << std::endl;

    _printDotRecursive(child, os);
  }
}

// use me for debugging :-)
template<typename T>
void Tree<T>::_printNodes() {
  for (unsigned i=0; i<m_vecNodes.size(); i++) { fprintf(stderr, "|%2u ", i); }; fprintf(stderr, "|\n");
  for (unsigned i=0; i<m_vecNodes.size(); i++) { fprintf(stderr, "+---"); }; fprintf(stderr, "+\n");
  for (unsigned i=0; i<m_vecNodes.size(); i++) {
    if (m_vecNodes[i] > 0) { fprintf(stderr, "|%s ", m_vecNodes[i]->label.c_str()); }
    else { fprintf(stderr, "| - "); }
  };
  fprintf(stderr, "|\n");
}

/*--------------------------------------*/
/*            parser logic              */
/*--------------------------------------*/

BOOST_FUSION_ADAPT_STRUCT(
  treeio::node,
  (treeio::children_vector, children)
  (string, label)
  (double, length)
)

namespace treeio {

  namespace qi = boost::spirit::qi;
  namespace ascii = boost::spirit::ascii;

  /** Represents the Newick format grammar.
    *
    * For details see http://evolution.genetics.washington.edu/phylip/newick_doc.html
    * All nodes must have labels.
    */
  template <typename Iterator>
  struct newick_grammar : qi::grammar<Iterator, node(), ascii::space_type>
  {
    newick_grammar() : newick_grammar::base_type(tree)
    {
      // use %= to assign the result of the parse to the string
      label %= quoted_label | unquoted_label;
      // unquoted labels may not contain formatting characters
      unquoted_label = qi::lexeme[*(qi::char_ - ':' - ',' - ';' - '(' - ')' - '[' - ']' - '\'')];
      // quoted labels can contain any characters
      quoted_label = '\'' >> qi::lexeme[*(qi::char_ - '\'')] >> '\'';

      // use %= to assign the result of the parse to the double
      branch_length %= ':' >> qi::double_;

      // Parser return values will be assigned to corresponding struct elements
      // in the order specified!
      subtree = -descendant_list >> -label >> -branch_length;

      // Descendant list is a vector of node, we just push back the
      // created nodes into the vector
      descendant_list = '(' >> subtree >> *(',' >> subtree ) >> ')';

      tree %= subtree >> ';';

      // prepare parsers for debugging
      BOOST_SPIRIT_DEBUG_NODE(label);
      BOOST_SPIRIT_DEBUG_NODE(branch_length);
      BOOST_SPIRIT_DEBUG_NODE(subtree);
      BOOST_SPIRIT_DEBUG_NODE(descendant_list);
      BOOST_SPIRIT_DEBUG_NODE(tree);
      // activate debugging output for parsers
      //qi::debug(tree);
      //qi::debug(subtree);
      //qi::debug(descendant_list);
      //qi::debug(label);
      //qi::debug(branch_length);
    }

    private:
      /* grammar rules, typed by element they create */
      qi::rule<Iterator, node(), ascii::space_type> tree, subtree;
      qi::rule<Iterator, vector<node>(), ascii::space_type> descendant_list;
      qi::rule<Iterator, string(), ascii::space_type> label;
      qi::rule<Iterator, string(), ascii::space_type> quoted_label;
      qi::rule<Iterator, string(), ascii::space_type> unquoted_label;
      qi::rule<Iterator, double(), ascii::space_type> branch_length;
  };

  template <typename Iterator>
  bool parse_newick(Iterator first, Iterator last, node &root_node)
  {
    newick_grammar<Iterator> grammar;

    bool result = qi::phrase_parse(first, last, grammar, ascii::space, root_node);
    if (first != last) return false;

    return result;
  }

  /***************************** readNewick *********************************/
  /** Reads a tree in Newick format */
  void readNewick (string newick_fn, node &tree) {
    ifstream f_newick;
    f_newick.open(newick_fn.c_str(), ios::in);
    readNewick(f_newick, tree);
    f_newick.close();
  }

  /***************************** readNewick *********************************/
  /** Reads a tree in Newick format */
  void readNewick (istream &input, node &tree) {
    string line;
    getline(input, line);
    if (parse_newick(line.begin(), line.end(), tree) == false) {
      throw runtime_error("Could not parse file in Newick format.");
    }
  }
}


// instantiate usable classes so they can be picked up by the linker
template class Tree<Clone>;
