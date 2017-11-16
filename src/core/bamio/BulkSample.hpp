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

  /**  */

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
  getCloneLabels ();

};

} // namespace bamio

#endif /* BULKSAMPLE_H */