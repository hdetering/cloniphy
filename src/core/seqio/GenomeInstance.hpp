#ifndef GENOMEINSTANCE_H
#define GENOMEINSTANCE_H

#include "AlleleSpecCopyNum.hpp"
#include "GenomeReference.hpp"
#include "ChromosomeInstance.hpp"

namespace seqio {

/** Represents the incarnation of a reference genome.
 *
 *  A GenomeInstance consists of a set of ChromosomeInstances,
 *  which in turn consist of SegmentCopies.
 */
struct GenomeInstance {
  /** Stores the ChromosomeInstances that make up this GenomeInstance. */
  std::vector<std::shared_ptr<ChromosomeInstance>> vec_chr;
  /** Stores for each ChromosomeInstance its real length. Used to randomly select a chromosome by length. */
  std::vector<unsigned long> vec_chr_len;
  /** ChromosomeInstances indexed by reference id. */
  std::map<std::string, std::vector<std::shared_ptr<ChromosomeInstance>>> map_id_chr;
  /** SegmentCopy interval maps indexed by chromosome id. */
  std::map<std::string, TSegMap> map_chr_seg;

  /** default c'tor */
  GenomeInstance();

  /** Create GenomeInstance from GenomeReference object
   *
   *  Create a diploid genome by initializing 2 ChromosomeInstances
   *  for each ChromosomeReference.
   */
  GenomeInstance(const GenomeReference& g_ref);

  /**
   * Create GenomeInstance as copy of existing object.
   * \param g_inst          Existing GenomeInstance to be copied.
   * \param out_vec_seg_mod Output param: mapping from original to new segment copy IDs.
   */
  GenomeInstance(
    const GenomeInstance& g_inst,
    std::vector<seg_mod_t>& out_vec_seg_mod
  );

  /** Adds a ChromosomeInstance to this GenomeInstance.
   *  \param sp_chr shared pointer to ChromosomeInstance
   *  \param id_chr reference ID under which to index ChromsomeInstance
   */
  void
  addChromosome(
    std::shared_ptr<ChromosomeInstance> sp_chr,
    std::string id_chr
  );

  /** Remove a ChromosomeInstance from this GenomeInstance.
   *  \param sp_chr       shared pointer to ChromosomeInstance
   *  \param id_chr       reference ID under which to index ChromsomeInstance
   *  \param idx_chr_cpy  index of ChromosomeInstance to delete
   */
  void
  deleteChromosomeInstance(
    std::shared_ptr<ChromosomeInstance> sp_chr,
    std::string id_chr,
    size_t idx_chr_cpy
  );

  /**
   * Perform Whole Genome Duplication.
   *
   * Duplication is carried out by duplicating all ChromosomeInstances.
   * \param out_vec_seg_mod Output parameter: vector of SegmentCopy modifications (id_new, id_old, start_old, end_old)
   */
  void
  duplicate(std::vector<seg_mod_t>& out_vec_seg_mod);

  /** Identify SequenceCopies overlapping a given locus. */
  std::vector<SegmentCopy>
  getSegmentCopiesAt (
    std::string id_chromosome,
    unsigned long ref_pos
  );

  /**DEPRECATED!
   * Map genomic regions to their absolute copy number states.
   *  \param out_map_cn_loci map storing for each copy number state a collection of genomic loci.
   */
  void
  getCopyNumberStates (
    std::map<unsigned, std::vector<std::shared_ptr<Locus>>>& out_map_cn_loci
  );

  /** Map genomic regions for each chromosome to their absolute copy number states.
   *  \param map_chr_seg map storing for genomic segment its copy number state.
   *  \param scale       factor by which each segment increases CN (default: 1).
   */
  void
  getCopyNumberStateByChr (
    std::map<std::string, boost::icl::interval_map<TCoord, double>>& map_chr_seg,
    const double scale = 1.0
  ) const;

  /** Map genomic regions for each chromosome to their allele-specific copy number states.
   *  \param map_chr_seg map storing for genomic segment its copy number state.
   *  \param scale       factor by which each segment increases CN (default: 1).
   */
  void
  getCopyNumberStateByChr (
    std::map<std::string, boost::icl::interval_map<TCoord, AlleleSpecCopyNum>>& map_chr_seg,
    const double scale = 1.0
  ) const;

  /** Build interval map of SegmentCopies for each chromosome.
   *  Used to answer questions like: 
   *    "Which SegmentCopies overlap interval 'chr1:10,000-10,500'?"
   *  \returns true on success, false on error
   */
   bool
   indexSegmentCopies ();
};

/** Print GenomeInstance substructure (for debugging purposes). */
std::ostream& operator<<(std::ostream& lhs, const GenomeInstance& gi);

} // namespace seqio

# endif // GENOMEINSTANCE_H