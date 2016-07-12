#ifndef BAMIO_H
#define BAMIO_H

#include "clone.hpp"
#include "random.hpp"
#include "vario.cpp"
#include <boost/format.hpp>
#include <memory>
#include <seqan/bam_io.h>
#include <string>

namespace bamio {

void mutateReads(
  std::string fn_fq_out,
  std::string fn_sam_out,
  std::string fn_sam_in,
  std::vector<vario::Variant> &variants,
  treeio::Tree<Clone> &tree,
  RandomNumberGenerator<> &rng
);

} // namespace bamio

#endif /* BAMIO_H */
