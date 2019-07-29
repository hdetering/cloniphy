#include "SegmentCopy.hpp"

using namespace std;

namespace seqio {

SegmentCopy::SegmentCopy ()
: id(boost::uuids::random_generator()()) {}

SegmentCopy::SegmentCopy (int i) {}

SegmentCopy::~SegmentCopy () {}

SegmentCopy::SegmentCopy (
  const TCoord start, 
  const TCoord end,
  const char allele
)
: id(boost::uuids::random_generator()()), 
  ref_start(start), 
  ref_end(end),
  gl_allele(allele)
{}

ostream& operator<<(ostream& lhs, const SegmentCopy& seg) {
  lhs << "SegmentCopy<uuid=" << seg.id << "> ";
  lhs << "[" << seg.ref_start << ", " << seg.ref_end << ")";
  lhs << endl;
  return lhs;
}

bool
operator==(const SegmentCopy& lhs, const SegmentCopy& rhs) {
  return lhs.id == rhs.id;
}

bool
operator<(const SegmentCopy& lhs, const SegmentCopy& rhs) {
  return lhs.id < rhs.id;
}

} // namespace seqio