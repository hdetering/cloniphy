#pragma once

#include "../seqio/types.hpp"
#include <string>

namespace vario {

/** 
 * CopyNumberVariants encapsulate CNV events with all of their properties 
 */
struct CopyNumberVariant
{
  unsigned    id;              /** unique identifier */
  bool        is_wgd;          /** true: CNV is Whole Genome Duplication */
  bool        is_deletion;     /** true: CNV is deletion event */
  bool        is_chr_wide;     /** true: event affects whole chromosome */
  bool        is_telomeric;    /** true: event coordinates include chromosome end */
  bool        is_first_arm;    /** true: first chr arm is affected; false: second arm. */
  bool        is_downstream;   /** true: event at 3' side of start_rel; false: at 5' end */
  double      len_rel;         /** length of affected region (fraction of chromsome length) */
  double      start_rel;       /** start position of event (fraction of chromosome length) */
  seqio::TCoord len_abs;       /** length of affected region (base pairs) */
  seqio::TCoord ref_pos_begin; /** start coordinate (in reference chr) */
  seqio::TCoord ref_pos_end;   /** end coordinate (in reference chr) */
  std::string ref_chr;         /** affected chromosome (reference ID) */
  std::string ref_allele;      /** affected allele ('A': maternal, 'B': paternal) */

  /** default c'tor */
  CopyNumberVariant();
};

} // namespace vario