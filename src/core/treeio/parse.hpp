#ifndef TREEIO_PARSE_H
#define TREEIO_PARSE_H

#include <iostream>
#include <string>
#include <vector>

namespace treeio {

/** Reads and writes tree files. */
namespace parse {

// forward declaration for later use
struct node;

// typedef to ease the writing
typedef std::vector<node> children_vector;

struct node
{
  std::string label;
  double length = 0;
  children_vector children;
};

void readNewick(std::string, node&);
void readNewick(std::istream&, node&);

} // namespace parse

} // namespace treeio

#endif // TREEIO_PARSE_H