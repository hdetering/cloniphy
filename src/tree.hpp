#ifndef TREE_H
#define TREE_H

#include "node.hpp"
#include <iostream>
#include <stdio.h>
#include <vector>

/** General tree representing a hierarchy of nodes. */
class Tree
{
  protected:
    int m_numTips;
    std::vector<Node *> m_vecNodes;
    Node *m_root;
    static void _printNewick(Node*, bool, std::ostream&);

  public:
    Tree(int);
    void generateRandomTopology(long);
    Node* getRoot();
    void printNodes();
    static void printNewick(Node*, std::ostream& = std::cout);
};

#endif /* TREE_H */
