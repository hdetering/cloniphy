#ifndef SEQUENCE_H
#define SEQUENCE_H

#include <string>

struct SequenceRecord
{
  std::string id;  /** Identifier */
  std::string seq; /** Actual sequence */
};

struct Mutation
{
  long absPos;  /** absolute bp position */
  int offset;   /** shifts the ancestral genotype to a new one */
};

#endif /* SEQUENCE_H */
