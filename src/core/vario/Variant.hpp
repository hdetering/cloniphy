#ifndef VARIANT_H
#define VARIANT_H

#include "../seqio/types.hpp"
#include <string>
#include <vector>

namespace vario {

/** Variants represent variable sites in nucleotide sequences */
struct Variant
{
  std::string id;     /** unique identifier */
  std::string chr;    /** reference chromosome id */
  short reg_copy;     /** affected copy of chr region (0: original, 1: first copy, ...) */
  seqio::TCoord pos;  /** reference basepair position */
  std::vector<std::string> alleles; /** observed alleles */
  int idx_mutation;   /** reference to mutation that gave rise to this variant */
  double rel_pos;     /** relative position in genome (use for sorting) */
  bool is_somatic;    /** is this Variant somatic or germline? (different output channels) */
  bool is_het;        /** true: Variant is heterozygous; false: homozygous (only applies if not is_somatic) */
  bool is_error;      /** true: Variant due to sequencing error (only applies to read count sim) */

  Variant();
  Variant(std::string id, std::string chr, unsigned long pos);
  ~Variant();

  bool operator< (const Variant&) const; /** make variants sortable */
  /** Sort variants by absolute position in genome. */
  static std::vector<Variant> sortByPosition(const std::vector<Variant>&);
  /** Sort variants, first lexicographically by CHR, then by position in CHR. */
  static std::vector<Variant> sortByPositionLex(const std::vector<Variant>&);
  /** Sort variants by absolute position in genome, taking chomosome copies into account. */
  //static std::vector<Variant> sortByPositionPoly(const std::vector<Variant>&);
  /** Sort variants by position in reference genome. */
  static std::vector<Variant> sortByPositionRef(const std::vector<Variant>&);
  /** Returns true if this variant is a SNV, false otherwise. */
  bool isSnv();
};

} // namespace vario

#endif // VARIANT_H