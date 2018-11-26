#ifndef SEQIO_TYPES_H
#define SEQIO_TYPES_H

#include <boost/uuid/uuid.hpp>
#include <boost/uuid/uuid_generators.hpp>
#include <boost/uuid/uuid_io.hpp>
#include <tuple>

namespace seqio {

/** Set of valid nucleotides. */
enum Nuc { 
  A, C, G, T, N 
};

/** Represents genomic coordinates. */
typedef
unsigned long 
TCoord;

/** Represents genomic regions (chromosome, start, end). */
typedef 
std::tuple<
  std::string, 
  TCoord, 
  TCoord
> 
TRegion;

/** Documents a modification introduced to a SegmentCopy.
 *  Structure: (new, old, start, end)
 *   - new: UUID of newly created SegmentCopy
 *   - old: UUID of existing SegmentCopy to be replaced
 *   - start: start coordinate within old SegmentCopy covered by the new one
 *   - end: end coordinate within old SegmentCopy covered by the new one
 */
typedef 
std::tuple<
  boost::uuids::uuid, 
  boost::uuids::uuid, 
  unsigned long, 
  unsigned long
>
seg_mod_t;


} // namespace seqio

#endif // SEQIO_TYPES_H