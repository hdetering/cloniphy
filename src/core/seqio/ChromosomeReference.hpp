#ifndef CHROMOSOMEREFERENCE_H
#define CHROMOSOMEREFERENCE_H

#include "SeqRecord.hpp"
#include "types.hpp"
#include <map>
#include <memory>

namespace seqio {

/** Represents the reference genome version of a chromosome.
 *
 *  Structurally a ChromosomeReference need not be a contiguous sequence,
 *  but may be composed of shorter elements (e.g., exons, targets).
 */
struct ChromosomeReference {
   /** chromosome identifier in reference genome. */
   std::string id;
   /** total length */
   unsigned long length;
   /** sequence records belonging to this chromosome, indexed by start position */
   std::map<TCoord, std::shared_ptr<SeqRecord>> map_start_rec;

   /** default c'tor */
   ChromosomeReference();

   /** Get nucleotide at position */
   char getNucAt(TCoord pos);
 };

} // namespace seqio

#endif // CHROMOSOMEREFERENCE_H