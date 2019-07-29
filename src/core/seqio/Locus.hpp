#ifndef LOCUS_H
#define LOCUS_H

#include "types.hpp"
#include <string>

namespace seqio {

/** Represents a genomic location */
struct Locus
{
  unsigned    idx_record; // index of genomic sequence record
  TCoord      start;      // global start position in sequence (0-based, inclusive)
  TCoord      end;        // global end position in sequence (0-based, exclusive)
  std::string id_ref;     // seq id in ref genome

  Locus();
  Locus(std::string id, TCoord start, TCoord end);
  ~Locus();
};

} // namespace seqio

#endif // LOCUS_H