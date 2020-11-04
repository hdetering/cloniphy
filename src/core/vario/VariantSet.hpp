#ifndef VARIANTSET_H
#define VARIANTSET_H

#include "Variant.hpp"
#include "../seqio.hpp"
#include <map>
#include <vector>

namespace vario {

/** VariantSets store a set of variants and summary statistics about them. */
struct VariantSet
{
  unsigned long num_variants = 0;
  std::vector<Variant> vec_variants; /** all variants that belong to the set */
  std::map<std::string, std::map<unsigned long, std::vector<Variant>>> map_chr2pos2var; /** variants stored by chromosome id */

  /** summary statistics */
  double mat_freqs[4][4] = {
    { 0.0, 0.0, 0.0, 0.0 },
    { 0.0, 0.0, 0.0, 0.0 },
    { 0.0, 0.0, 0.0, 0.0 },
    { 0.0, 0.0, 0.0, 0.0 }
  }; /** nucleotide substitution frequencies */

  VariantSet();
  VariantSet(std::vector<Variant> variants);
  ~VariantSet();

  /** Compound assignment VariantSets: add Variants of rhs to this VariantSet.
   *  (Does not need to be a member, but often is, to modify the private members)
   */
  VariantSet& operator+=(const VariantSet& rhs);
  /** Add two VariantSets: returns the union set.
   *  Friends defined inside class body are inline and are hidden from non-ADL lookup.
   *  Passing lhs by value helps optimize chained a+b+c. (Otherwise, both parameters
   *    may be const references)
   */
  friend VariantSet operator+(VariantSet lhs, const VariantSet& rhs) {
    lhs += rhs; // reuse compound assignment
    return lhs; // return the result by value (uses move constructor)
  }

  /** Index variants by chromosome and position. */
  long indexVariants();
  long calculateSumstats();
};

} // namespace vario

#endif // VARIANTSET_H