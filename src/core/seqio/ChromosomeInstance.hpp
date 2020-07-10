#ifndef CHROMOSOMEINSTANCE_H
#define CHROMOSOMEINSTANCE_H

#include "ChromosomeReference.hpp"
#include "SegmentCopy.hpp"
#include "types.hpp"
#include <boost/icl/interval.hpp>
#include <boost/icl/interval_map.hpp>
#include <boost/uuid/uuid.hpp>
#include <boost/uuid/uuid_generators.hpp>
#include <boost/uuid/uuid_io.hpp>
#include <list>

namespace seqio {

/** Stores SegmentCopies overlapping genomic intervals. */
typedef 
boost::icl::interval_map<
  TCoord, 
  std::set<SegmentCopy>
> 
TSegMap;

/** Stores SegmentCopies as interval set */
typedef 
boost::icl::interval_set<TCoord> 
TSegSet;

/** Represents a physical representation of a chromosome.
 *
 *  Structurally a ChromosomeInstance is a chain of SegmentCopies,
 *  which represent regions of the reference sequence.
 */
struct ChromosomeInstance {
  /** Universally-unique identifier to globally refer to the object. */
  boost::uuids::uuid id;
  /** Identifier of chromosome in reference genome */
  std::string id_ref;
  /** Real length (in bp) of ChromosomeInstance */
  unsigned long length;
  /** SegmentCopies that are associated with this ChromosomeInstance */
  std::list<SegmentCopy> lst_segments;
  
  /** default c'tor */
  ChromosomeInstance();

  /** Generate ChromosomeInstance from reference chromosome.
   *  A single SegmentCopy will be created comprising the whole chromosome.
   *  \param chr_ref    Reference chromosome of which to create an instance.
   *  \param gl_allele  Germline source allele, will be assigned to segment copies.
   */
  ChromosomeInstance (
    const ChromosomeReference chr_ref,
    const char gl_allele
  );

  /** Identify SegmentCopies overlapping a given locus. */
  std::vector<SegmentCopy> 
  getSegmentCopiesAt (
    TCoord ref_pos
  );

  /** Amplify a region of this ChromosomeInstance by creating new SegmentCopies
   *  \param start_abs    Output param: absolute start coord in reference (bp)
   *  \param end_abs      Output param: absolute end coord in reference (bp)
   *  \param len_abs      Output param: absolute length of event (bp)
   *  \param start_rel    Relative start coordinate of deletion (fraction of chromosome length).
   *  \param len_rel      Relative length of region to delete (fraction of chromosome length).
   *  \param is_forward   true: Deletion affects region towards 3' end from start_rel; false: towards 5' end.
   *  \param is_telomeric true: Deletion includes chromosome end (facilitates end coordinate calculation)
   *  \param is_first_arm true: First chrom arm is affected; false: second arm affected
   */
  std::vector<seg_mod_t> 
  amplifyRegion (
    seqio::TCoord& start_abs,
    seqio::TCoord& end_abs,
    seqio::TCoord& len_abs,
    const double start_rel,
    const double len_rel,
    const bool is_forward,
    const bool is_telomeric,
    const bool is_first_arm
  );

  /** Delete a region of this ChromosomeInstance by creating new SegmentCopies.
   *  \param start_abs    Output param: absolute start coord in reference (bp)
   *  \param end_abs      Output param: absolute end coord in reference (bp)
   *  \param len_abs      Output param: absolute length of event (bp)
   *  \param start_rel    Relative start coordinate of deletion (fraction of chromosome length).
   *  \param len_rel      Relative length of region to delete (fraction of chromosome length).
   *  \param is_forward   true: Deletion affects region towards 3' end from start_rel; false: towards 5' end.
   *  \param is_telomeric true: Deletion includes chromosome end (facilitates end coordinate calculation)
   *  \param is_first_arm true: First chrom arm is affected; false: second arm affected
   */
  std::vector<seg_mod_t> 
  deleteRegion (
    seqio::TCoord& start_abs,
    seqio::TCoord& end_abs,
    seqio::TCoord& len_abs,
    const double start_rel,
    const double len_rel,
    const bool is_forward,
    const bool is_telomeric,
    const bool is_first_arm
  );

  /** Create a copy of an existing ChromosomeInstance.
   *
   *  Copies each SegmentCopy comprising the existing ChromosomeInstance.
   *  \param ci_old Existing ChromosomeInstance to be copied
   *  \param out_seg_mods Output parameter: Tuples of modifications to SegmentCopies (id_new, id_old, start_old, end_old)
   */
  void 
  copy (
    std::shared_ptr<ChromosomeInstance> ci_old, 
    std::vector<seg_mod_t>& out_seg_mods
  );

  /** Build interval map of SegmentCopies.
   *  Used to answer questions like: 
   *    "Which SegmentCopies overlap interval 'chr1:10,000-10,500'?"
   *  \param imap_segments output param: interval map of SegmentCopies
   *  \returns true on success, false on error
   */
   bool
   indexSegmentCopies (
     TSegMap& imap_segments
  ) const;

}; // struct ChromosomeInstance

/** Print ChromosomeInstance substructure (for debugging purposes). */
std::ostream& operator<<(std::ostream& lhs, const ChromosomeInstance& ci);

} // namespace seqio

#endif // CHROMOSOMEINSTANCE_H