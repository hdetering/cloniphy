#ifndef SEQRECORD_H
#define SEQRECORD_H

#include <string>

namespace seqio {

struct SeqRecord
{
  std::string id;          /** identifier */
  std::string description; /** sequence description (everything after first space in ID line) */
  std::string seq;         /** actual sequence */
  std::string id_ref;      /** identifier in reference genome (ploidy) */
  short chr_copy;          /** chromosome copy (0 for haploid) */
  SeqRecord(const std::string, const std::string, const std::string&);
  ~SeqRecord();
};

} // namespace seqio

#endif // SEQRECORD_H