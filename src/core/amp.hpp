#ifndef AMP_H
#define AMP_H

#include "random.hpp"
#include <vector>

namespace amp {

//! Amplicons represent physical copies of a stretch of DNA.

struct Amplicon {
  //! Copy of a chromosome that Amplicon represents (diploid: {0,1})
  unsigned short allele;
  //! Amplicon physical length
  unsigned length;
};

//! Simulate Multiple Displacement Amplification (MDA).
/*!
  \param cvg output cvg vector of amplicon coverage values (index is position in sequence)
  \param seq_len length of primary sequence
  \param amplicon_size maximum size of generated amplicons
  \param fold factor by which the sequence should be amplified (determines when to stop)
  \param rng random number generator
 */
void simulateMda(
  std::vector<unsigned> &cvg,
  unsigned seq_len,
  double amplicon_size_mean,
  double amplicon_size_sd,
  int fold,
  RandomNumberGenerator &rng
);

//! Simulate Multiple Displacement Amplification (MDA) of a single DNA molecule.
/*!
  \param cvg output cvg vector of amplicon coverage values (index is position in sequence)
  \param seq_len length of primary sequence
  \param amplicon_size maximum size of generated amplicons
  \param fold factor by which the sequence should be amplified (determines when to stop)
  \param rng random number generator
 */
void simulateMdaProcessHaploid(
  std::vector<unsigned> &cvg,
  unsigned seq_len,
  double amplicon_size_mean,
  double amplicon_size_sd,
  int fold,
  RandomNumberGenerator &rng
);

//! Simulate Multiple Displacement Amplification (MDA) of a diploid chromosome.
/*!
  \param cvg output vector of amplicon coverage vectors (one for each allele, index is position in sequence)
  \param seq_len length of primary sequence
  \param amplicon_size maximum size of generated amplicons
  \param fold factor by which the sequence should be amplified (determines when to stop)
  \param rng random number generator
 */
void simulateMdaProcessDiploid(
  std::vector<std::vector<unsigned>> &cvg,
  unsigned long long seq_len,
  double amplicon_size_mean,
  double amplicon_size_sd,
  int fold,
  RandomNumberGenerator &rng
);


} /* namespace amp */

#endif /* AMP_H */
