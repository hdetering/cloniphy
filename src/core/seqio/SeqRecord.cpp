#include "SeqRecord.hpp"

using namespace std;

namespace seqio {

SeqRecord::SeqRecord(const string id, const string desc, const string& seq)
  : id(id), description(desc), seq(seq), id_ref(id), chr_copy(0) {}
SeqRecord::~SeqRecord() {}

} // namespace seqio