#ifndef BEDFILE_H
#define BEDFILE_H

#include "types.hpp"
using seqio::TCoord;
#include "Locus.hpp"
using seqio::Locus;
#include "../stringio.hpp"
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

namespace seqio {

/** 
 * Handles input, output and conversion of BED files.
 * 
 * A BED file is expected to have three required columns:
 *   1. seq_id : Identifier for a sequence (e.g. chromosome)
 *   2. start  : Start coordinate (0-based, inclusive)
 *   3. end    : End coordinate (0-based, exclusive)
 * 
 * Additional columns are collected in a feature collection.
 */
struct BedFile {

/** Number of records in this BED file. */
size_t m_num_recs;
/** Number of additional features (apart from seqid, start, end). */
size_t m_num_feat;
/** Sequence IDs. */
std::vector<std::string> m_vec_seqid;
/** Start coordinates (0-based, inclusive). */
std::vector<unsigned long> m_vec_start;
/** End coordinates (0-based, exclusive). */
std::vector<unsigned long> m_vec_end;
/** Additional columns' values. */
std::vector<std::vector<std::string>> m_feat_val;
/** Additional columns' names. */
std::vector<std::string> m_feat_names;

/** default c'tor */
BedFile();
/** Initialize from existing BED file. */
BedFile(const std::string& filename);
/** Initialize from existing BED file. */
BedFile (
  const std::string& filename, 
  const std::vector<std::string>& feat_names
);

/** Get records from BED file. */
bool parseBed(const std::string& filename);
/** Get records from BED file. */
bool parseBed(std::ifstream& input);

/** Return genomic records as list of Locus objects. */
void getLoci(std::vector<Locus>& out_loci);

};

} /* namespace seqio */

#endif /* BEDFILE_H */