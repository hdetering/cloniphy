#ifndef CLONE_H
#define CLONE_H

#include "node.hpp"

struct Clone: public Node
{
  float freq;
  bool is_healthy;
};

#endif /* CLONE_H */
