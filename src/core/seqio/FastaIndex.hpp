#ifndef FASTAINDEX_H
#define FASTAINDEX_H

#include "types.hpp"
using seqio::TCoord;
#include "../stringio.hpp"
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

namespace seqio {

/**
 * Reads and converts FASTA index files (as created by `samtools faidx`).  
 */
struct FastaIndex {

/** Number of sequences. */
size_t m_num_seqs;
/** List of sequence names. */
std::vector<std::string> m_vec_name;
/** List of sequence lengths. */
std::vector<unsigned long> m_vec_len;
/** Byte offset of sequence start in FASTA file. */
std::vector<unsigned long> m_vec_offset;
/** Sequence line length. */
std::vector<unsigned long> m_vec_linebases;
/** Sequence line length, including newline character. */
std::vector<unsigned long> m_vec_linewidth;

/** default c'tor */
FastaIndex();
/** Initialize from existing .fai file. */
FastaIndex(const std::string& fn_fai);

/** Read records from .fai file stream. */
bool parseFai(std::ifstream& input);

};

} /* namespace seqio */

#endif /* FASTAINDEX_H */