#ifndef BULKSAMPLE_H
#define BULKSAMPLE_H

#include "../bamio.hpp"

namespace bamio {

/**
 * Encapsulates information about a bulk mixture of clones. 
 */
struct BulkSample {
  
  /* PROPERTIES */

  /** sample ID */
  std::string id;

  /** clone weights (cell fraction) */
  std::map<std::string, double> m_clone_weight;

  /** absolute genome length (all segment copies considered) */
  seqio::TCoord genome_len_abs;

  /** Allele-specific copy number state by chromosome, coordinates */
  std::map<
    std::string,     // chromosome id
    boost::icl::interval_map<
      seqio::TCoord, // start and end coordinates 
      seqio::AlleleSpecCopyNum // copy no. for maternal and paternal allele
    >
  > m_chr_cn;

  /* METHODS */

  /** std c'tor */
  BulkSample ();

  /* Initializes a new BulkSample object from id and clone weights.
   * 
   * \param id_sample     Unique identifier for sample. 
   * \param clone_weight  Clone weights, represent cancer cell fraction (CCF).
   */
  BulkSample ( 
    const std::string id_sample,
    const std::map<std::string, double> clone_weight
  );

  /**
   * Returns clone labels of clones in this sample.
   * 
   * NOTE: Clones having cellular frequency of 0.0 will also be returned.
   * 
   * \returns  Clone labels belonging to sample.
   */
  std::vector<std::string>
  getCloneLabels () const;

  /** 
   * Return copy number-adjusted expected sequencing coverage for locus.
   * 
   * expected coverage for affected segment S_i is calculated as follows:
   *   c(i) = C_bp * n(i) * l(i)
   * with
   *   C_bp  := expected coverage per base pair
   *   n(i)  := combined copy number for genomic segment S_i
   *   l(i)  := length of genomic segment S_i
   * 
   * and
   *   C_bp = C * G_ref / G_abs
   * with
   *   C     := global target coverage (relative to ref. genome)
   *   G_ref := reference genome length
   *   G_abs := absolute genome length = \sum_i ( n(i) * l(i) )
   * 
   * \param chr         affected chromosome
   * \param pos         affected base pair position
   * \param cvg_per_bp  expected coverage per base pair
   */
  double 
  getExpectedCoverageAt (
    std::string chr,
    seqio::TCoord pos,
    double cvg_per_bp
  );

  /**
   * Return total copy number state for given genomic position.
   * 
   * Total copy number: sum of maternal and paternal allele copy number.
   * 
   * \param chr  Chromosome id.
   * \param pos  Reference base pair position in chromosome.
   * \param out_cn_tot   Out param: total copy number at given locus.
   * \param out_seq_len  Out param: length of genomic segment harboring given locus.
   * \returns    true on success, false on error
   */
  bool
  getTotalCopyNumberAt (
    const std::string chr,
    const seqio::TCoord pos,
    double out_cn_tot,
    seqio::TCoord out_seg_len
  ) const;

};

} // namespace bamio

#endif /* BULKSAMPLE_H */