#ifndef NODE_H
#define NODE_H
#include <vector>

/** Generic tree node */
struct Node
{
  int label;
  Node *parent;

  Node();
  virtual ~Node();
  virtual float distanceToParent() = 0;
  virtual bool isLeaf() = 0;
  bool isRoot();
};

#endif /* NODE_H */
