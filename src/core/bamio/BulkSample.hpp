#ifndef BULKSAMPLE_H
#define BULKSAMPLE_H

#include "../bamio.hpp"
#include "../seqio/AlleleSpecCopyNum.hpp"

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

  /** List of all genomic segments contained in this BulkSample.
   * Initialized by BulkSampleGenerator::calculateBulkCopyNumber.
   * Used to select locations of sequencing errors during read count simulation.
   * Using seqio::Locus here because it knows which chromosome it belongs to.
   */
  std::vector<seqio::Locus> m_vec_loc;
  /** Weights for SegmentCopies in m_vec_seg.
   * Initialized by BulkSampleGenerator::calculateBulkCopyNumber.
   * Used to select locations of sequencing errors during read count simulation.
   * Weight is calculated by segment length * mixed copy number.
   */
  std::vector<double> m_vec_loc_weight;

  /** Allele counts of SNVs indexed by clone, SNV id. */
  std::map<std::string, std::map<int, vario::VariantAlleleCount>> m_map_clone_snv_vac;
  /** Allele frequencies of SNVs indexed by SNV id. */
  std::map<int, double> m_map_snv_vaf;

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
   * Get total copy number, depending on clonal mixture of sample.
   * 
   * \param chr   affected chromosome
   * \param pos   affected base pair position
   * \returns     copy number for sample at locus chr:pos
   */
  // double
  // getTotalCopyNumberAt (
  //   const std::string chr,
  //   const seqio::TCoord pos
  // ) const;

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
    const std::string chr,
    const seqio::TCoord pos,
    const double cvg_per_bp
  ) const;

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
    double& out_cn_tot,
    seqio::TCoord& out_seg_len
  ) const;

  /** 
   * Calculate allele counts at variant positions for all clones.
   * Initializes m_map_clone_snv_vac, m_map_snv_vaf.
   *
   * \param map_clone_ccf      Cancer cell fraction (CCF) for each clone in the sample.
   * \param var_store          Variant store keeping track of SNV -> segment copy mappings.
   * \param map_clone_chr_seg  Genomic segment copies for each chromosome and clone.
   * \returns                  True on success, false on error.
   */
  bool
  initAlleleCounts (
    const std::map<std::string, double> map_clone_ccf,
    const vario::VariantStore& var_store,
    const std::map<std::string, std::map<std::string, seqio::TSegMap>> map_clone_chr_seg
  );
};

} // namespace bamio

#endif /* BULKSAMPLE_H */