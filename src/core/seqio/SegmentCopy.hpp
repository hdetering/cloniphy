#ifndef SEGMENTCOPY_H
#define SEGMENTCOPY_H

#include "types.hpp"
#include <boost/uuid/uuid.hpp>
#include <boost/uuid/uuid_generators.hpp>
#include <boost/uuid/uuid_io.hpp>

namespace seqio {

/** Represents a copy of a reference sequence segment.
 *
 *  Inititally, each allele (n=2 in the normal, diploid case) is a copy of
 *  the reference sequence. Further segment copies arise from CNV events.
 */
struct SegmentCopy {
  /** Universally-unique identifier to globally refer to the object. */
  boost::uuids::uuid id;
  /** Start position (0-based, inclusive) in reference chromosome */
  TCoord ref_start;
  /** End position (0-based, exclusive) in reference chromosome */
  TCoord ref_end;
  /** Germline source allele (maternal: A, paternal: B). */
  char gl_allele;

  /** default c'tor */
  SegmentCopy();
  /** dummy c'tor to avoid UUID init (to avoid performance overhead) */
  SegmentCopy(int);
  /** default d'tor */
  ~SegmentCopy();

  /** Create SegmentCopy for given coordinates.
   *  \param start      Reference bp position where segment copy begins.
   *  \param end        Reference bp position where segment copy ends.
   *  \param gl_allele  Germline source allele ('A','B') from which segment copy originates.
   */
  SegmentCopy (
    const TCoord start, 
    const TCoord end,
    const char gl_allele
  );
};

/** Compare two segment copies (equality).
 *  If IDs match, segment copies are identical.
 */
bool operator==(const SegmentCopy& lhs, const SegmentCopy& rhs);

/** Compare two segment copies (order).
 *  Order by UUID (necessary to avoid collisions in boost::icl::set).
 */
bool operator<(const SegmentCopy& lhs, const SegmentCopy& rhs);

/** Print SegmentCopy info (for debugging purposes). */
std::ostream& operator<<(std::ostream& lhs, const SegmentCopy& seg);

} // namespace seqio

#endif // SEGMENTCOPY_H