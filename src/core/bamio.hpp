#ifndef BAMIO_H
#define BAMIO_H

#include "clone.hpp"
#include "random.hpp"
//#include "vario.cpp"
#include "stringio.hpp"
#include <boost/format.hpp>
#include <memory>
#include <seqan/bam_io.h>
#include <string>

/**
 * Methods to support the input/output of SAM/BAM files.
 */
namespace bamio {

/** Takes an existing SAM/BAM file and adds germline mutations to reads. */
void mutateReads(
  const std::string fn_sam_out,
  const std::string fn_sam_in,
  const vario::VariantSet &variants,
  const short ploidy,
  RandomNumberGenerator<> &rng
);

/** Takes an existing SAM/BAM file and adds subclonal mutations to reads.
 *  NOTE: weights are to be provided for clone nodes (root node receives 1-sum)
 */
void mutateReads(
  std::string fn_fq_out,
  std::string fn_sam_out,
  std::string fn_sam_in,
  vario::VariantSet &variants,
  treeio::Tree<Clone> &tree,
  std::vector<double> weights,
  std::string id_sample,
  const short ploidy,
  RandomNumberGenerator<> &rng,
  bool do_write_fastq = false
);

/** Add a read group to BAM file header for each clone. */
void addCloneReadGroups(
  seqan::BamHeader &header,
  const std::string id_sample,
  const std::vector<std::string> &vec_lbl
);

} // namespace bamio

#endif /* BAMIO_H */
