#ifndef AMP_H
#define AMP_H

#include "random.hpp"
#include <vector>

namespace amp {

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
  RandomNumberGenerator<> rng
);

//! Simulate Multiple Displacement Amplification (MDA).
/*!
  \param cvg output cvg vector of amplicon coverage values (index is position in sequence)
  \param seq_len length of primary sequence
  \param amplicon_size maximum size of generated amplicons
  \param fold factor by which the sequence should be amplified (determines when to stop)
  \param rng random number generator
 */
void simulateMdaProcess(
  std::vector<unsigned> &cvg,
  unsigned seq_len,
  double amplicon_size_mean,
  double amplicon_size_sd,
  int fold,
  RandomNumberGenerator<> rng
);

} /* namespace amp */

#endif /* AMP_H */
