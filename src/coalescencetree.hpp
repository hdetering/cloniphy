#include "node.hpp"

/** Binary tree representing the genealogy of samples. */
class CoalescenceTree
{
  // if you don't declare private members in the header -> seqfault on destruction!
  std::vector<Node> m_vecNodes;
  Node *m_root;

  public:
    CoalescenceTree(int);
    void randomize();
    Node getRoot();

};
