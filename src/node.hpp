#ifndef NODE_H
#define NODE_H

/** Generic tree node */
struct Node
{
  int label;
  Node *parent;
  Node *left;
  Node *right;

  Node();
  bool isLeaf();
  bool isRoot();
};

#endif /* NODE_H */
