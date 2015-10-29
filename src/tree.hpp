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

  public:
    virtual Node* getRoot() = 0;

};

#endif /* TREE_H */
