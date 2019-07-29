#ifndef ALLELESPECCOPYNUM_H
#define ALLELESPECCOPYNUM_H

namespace seqio {

/** Utility class to encapsulate allele-specific copy number state as two double values. */
struct AlleleSpecCopyNum {
  /** Copy number for allele A. */
  double count_A;
  /** Copy number for allele B. */
  double count_B;

  /** default c'tor */
  AlleleSpecCopyNum() : count_A(0.0), count_B(0.0) {}

  /** Add two copy number states. */
  AlleleSpecCopyNum& operator+=(const AlleleSpecCopyNum& rhs) {
    this->count_A += rhs.count_A;
    this->count_B += rhs.count_B;
    return *this;
  }
};
/** Compare two copy number states. */
bool operator==(const AlleleSpecCopyNum& lhs, const AlleleSpecCopyNum& rhs);

} // namespace seqio

#endif // ALLELESPECCOPYNUM_H