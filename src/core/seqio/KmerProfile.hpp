#ifndef KMERPROFILE_H
#define KMERPROFILE_H

#include "../stringio.hpp"
#include <map>
#include <string>
#include <vector>

namespace seqio {

/**
 * Encapsulates trinucleotide profile for a reference sequence.
 * Each trinucleotide (AAA, AAC, ..., TTG, TTT) is assigned a probabilty/weight.
 */
struct KmerProfile {

/** Length of kmers. */
short kmer_length = 0;
/** Number of kmers in index. */
unsigned num_kmers = 0;
/** Stores list of kmers. */
std::vector<std::string> 
  m_vec_kmers;
/** Stores list of kmer weights. Position corresponds to m_vec_kmers. */
std::vector<double>
  m_vec_weight;
/** Stores weights for prefixes of length k-1. */
std::map<std::string, KmerProfile> m_idx_pfx;
/** Flag indicating that prefices have been indexed. */
bool has_idx_pfx;

/** default c'tor */
KmerProfile();
/** Initialize Profile from an existing CSV file. */
KmerProfile(std::string fn_csv);

/** Read trinucleotide from a CSV file. */
bool
readCSV (
  const std::string filename
);

/** Calculate sub-indices for k-1 prefixes. */
bool indexPrefixes ();

}; // struct KmerProfile

} // namespace seqio

#endif // KMERPROFILE_H