#ifndef NODE_H
#define NODE_H
#include <vector>

/** Generic tree node */
struct Node
{
  int label;
  Node *parent;
  std::vector<Node*> m_vecChildren;

  Node();
  bool isLeaf();
  bool isRoot();
};

#endif /* NODE_H */
