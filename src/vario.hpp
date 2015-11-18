#ifndef VARIO_H
#define VARIO_H

#include "seqio.hpp"

#include <iostream>
#include <vector>

/** Mutations specifiy a modification of a sequence.
 *  TODO: Do I really belong here?
 */
struct Mutation
{
  unsigned id;          /** unique identifier */
  unsigned long absPos; /** absolute bp position */
  short offset;         /** shifts the ancestral genotype to a new one */
  short copy;           /** which chromosome copy (0:maternal, 1:paternal) */

  bool operator< (const Mutation&) const; /** make mutations sortable */
  static std::vector<Mutation> sortByPosition(const std::vector<Mutation>&);
};

class VarIO {
  public:
    /** Generate VCF output from a reference sequence and a set of mutations.  */
    static void writeVcf(const std::vector<SeqRecord>&, const std::vector<Mutation>&, const std::vector<std::vector<short> > &, std::ostream&);
};

#endif /* VARIO_H */
