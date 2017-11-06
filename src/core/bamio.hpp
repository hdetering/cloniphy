#ifndef BAMIO_H
#define BAMIO_H

#include "clone.hpp"
#include "random.hpp"
#include "seqio.hpp"
#include "vario.hpp"
#include "stringio.hpp"
#include <boost/filesystem.hpp> // absolute(),
#include <boost/format.hpp>
#include <memory>
#include <seqan/bam_io.h>
#include <string>

using vario::Variant;

/**
 * Methods to support the input/output of SAM/BAM files
 */
namespace bamio {

/** encapsulates sequence simulation with ART */
class ArtWrapper {
public:
  std::string bin_path;
  unsigned read_len;
  unsigned frag_len_mean;
  unsigned frag_len_sd;
  float fold_cvg;
  /** number of read (pairs) to sample, NOTE: overrides fold_cvg! */
  unsigned long num_reads;
  std::string out_pfx;
  bool out_sam;
  bool out_aln;
  std::string seq_sys;
  /** path to reference FASTA file */
  std::string fn_ref_fa;
  /** should FASTQ output file be kept? */
  bool do_keep_fq;

  ArtWrapper(std::string bin_path);
  int run(std::string out_pfx);
};

/** Takes an existing SAM/BAM file and adds germline mutations to reads. */
void mutateReads(
  const std::string fn_sam_out,
  const std::string fn_sam_in,
  const vario::VariantSet &variants,
  const short ploidy,
  RandomNumberGenerator<> &rng
);

/** DEPRECATED Takes an existing SAM/BAM file and adds subclonal mutations to reads.
 *  NOTE: weights are to be provided for clone nodes (root node receives 1-sum)
 */
void 
mutateReads (
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

/** Takes an existing SAM/BAM file and adds subclonal mutations to reads.
 *  NOTE: weights are to be provided for clone nodes (root node receives 1-sum)
 */
void 
mutateReads (
  boost::filesystem::path fn_fq_out,
  boost::filesystem::path fn_sam_out,
  boost::filesystem::path fn_sam_in,
  vario::VariantSet &variants,
  std::vector<std::shared_ptr<Clone>> vec_clones,
  std::map<std::string, std::vector<bool>> mm,
  std::vector<double> weights,
  std::string id_sample,
  const short ploidy,
  RandomNumberGenerator<> &rng,
  bool do_write_fastq = false
);

/** Add a read group to BAM file header for each clone. */
void 
addCloneReadGroups (
  seqan::BamHeader &header,
  const std::string id_sample,
  const std::vector<std::string> &vec_lbl
);

} // namespace bamio

#endif /* BAMIO_H */
