#include "Locus.hpp"

using namespace std;

namespace seqio {
  
Locus::Locus() {}
Locus::Locus(string id, ulong start, ulong end)
: id_ref(id), start(start), end(end), idx_record(0) {}
Locus::~Locus() {}

} // namespace seqio