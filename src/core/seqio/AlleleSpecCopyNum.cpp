#include "AlleleSpecCopyNum.hpp"

namespace seqio {

bool 
operator==(const AlleleSpecCopyNum& lhs, const AlleleSpecCopyNum& rhs) {
  return lhs.count_A == rhs.count_A && lhs.count_B == rhs.count_B;
}

} // namespace seqio