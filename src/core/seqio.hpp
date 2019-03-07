#ifndef SEQIO_H
#define SEQIO_H

#include "seqio/KmerProfile.hpp"
#include "seqio/Locus.hpp"
#include "seqio/SegmentCopy.hpp"
#include "seqio/SeqRecord.hpp"
#include "seqio/types.hpp"
#include "random.hpp"
#include "stringio.hpp"
#include <algorithm>
#include <boost/circular_buffer.hpp>
#include <cassert>
#include <fstream>
#include <functional>
#include <iostream>
#include <map>
#include <memory> // unique_ptr, shared_ptr, weak_ptr
#include <sstream>
#include <string>
#include <vector>

/** Handles sequence files (read, write, index). */
namespace seqio {

/**
 * Returns a canonical index (corresponding to ::Nuc) for a nucleotide.
 *
 * \returns index position of nucleotide, -1 for unknown characters
 */
inline short nuc2idx (char nuc) {
  switch (nuc) {
    case 'A': case 'a': return 0;
    case 'C': case 'c': return 1;
    case 'G': case 'g': return 2;
    case 'T': case 't': return 3;
  }
  return -1; // default for unknown chars
}

/**
 * Map nucleotide char to a canonical index (corresponding to ::Nuc).
 *
 * \returns Nucleotide char for defined indices, `?` otherwise
 */
inline char idx2nuc (short idx) {
  switch (idx) {
    case 0: return 'A';
    case 1: return 'C';
    case 2: return 'G';
    case 3: return 'T';
    case 4: return 'N';
  }
  return '?'; // default for unknown chars
}

/**
 * Get reverse complement of a nucleotide.
 */
inline char rev_comp (char nuc) {
   switch (nuc) {
     case 'A': case 'a': return 'T';
     case 'C': case 'c': return 'G';
     case 'G': case 'g': return 'C';
     case 'T': case 't': return 'A';
   }
   assert(false);
   return '?'; // default for unknown chars
}

/**
 * Get reverse complement of a DNA sequence.
 */
inline std::string rev_comp (std::string dna) {
  std::transform(dna.begin(), dna.end(), dna.begin(),
                 [](char c) { return rev_comp(c); });
  return dna;
}

/** 
 * Read sequences from FASTA file. 
 * \param records   output parameter; list of records.
 * \param filename  filename to read records from.
 *
void readFasta (
  std::vector<std::shared_ptr<SeqRecord>>& records,
  const char* filename
);*/
typedef std::map<std::string, std::vector<Locus>> TLocusMap;
/** 
 * Read sequences from FASTA file. 
 * \param records        output parameter; list of records.
 * \param filename       filename to read records from.
 * \param use_whitelist  only import sequences specified in map_seq_loci.
 * \param map_seq_loci   sequence segments to import (whitelist).
 */
void readFasta (
  std::vector<std::shared_ptr<SeqRecord>>& records,
  const char* filename,
  const bool use_whitelist = false,
  const TLocusMap& map_seq_loci = TLocusMap()
);
/** 
 * Reads sequences from istream. 
 * \param records      output parameter; list of records.
 * \param inputstream  filename to read records from.
 *
void readFasta (
  std::vector<std::shared_ptr<SeqRecord>>& records,
  std::istream& inputstream
);*/
/** 
 * Reads sequences from istream. 
 * \param records        output parameter; list of records.
 * \param input_stream   input stream to read records from.
 * \param use_whitelist  only import sequences specified in map_seq_loci.
 * \param map_seq_loci   sequence segments to import (whitelist).
 */
void readFasta (
  std::vector<std::shared_ptr<SeqRecord>>& records,
  std::istream& inputstream,
  const bool use_whitelist = false,
  const TLocusMap& map_seq_loci = TLocusMap()
);

/** Writes sequences to file. */
int writeFasta(const std::vector<std::shared_ptr<SeqRecord>>&, const std::string fn, int len_line = 60);
/** Writes sequences to ostream. */
int writeFasta(const std::vector<std::shared_ptr<SeqRecord>>&, std::ostream& os, int len_line = 60);
/** Generate an index for a FASTA file containing multiple sequences */
void indexFasta(const char*);

/** Generate random DNA sequence with defined nucleotide frequencies. 
 *  \param seq        output parameter; receives the sequence that is generated.
 *  \param total_len  length of sequence to generate.
 *  \param nuc_freqs  nucleotide frequencies (A,C,G,T).
 *  \param rng        random number generator.
 *  \returns          length of generated sequence.
 */
unsigned long 
generateRandomDnaSeq (
  std::string &seq,
  const unsigned long total_len,
  const std::vector<double> nuc_freqs,
  RandomNumberGenerator &rng);

/** Generate random DNA sequence with defined kmer frequencies. 
 *  \param seq        output parameter; receives the sequence that is generated.
 *  \param total_len  length of sequence to generate.
 *  \param kmer_prof  Kmer frequency profile.
 *  \param rng        random number generator.
 *  \returns          length of generated sequence.
 */
unsigned long 
generateRandomDnaSeq (
  std::string &seq,
  const unsigned long total_len,
  const KmerProfile& kmer_prof,
  RandomNumberGenerator &rng);

/** Simulate allelic dropout events, masking parts of genome as 'N's. */
void simulateADO_old(const std::string, const float, const int, std::function<double()>&);
/** Simulate allelic dropout events, masking parts of genome as 'N's. */
void simulateADO(const std::string, const unsigned, const float, const int, std::function<double()>&);
/** Convert a nucleotide char into Nuc */
Nuc charToNuc(const char);
/** Convert Nuc into a nucleotide char */
char nucToChar(const Nuc);
/** Convert Nuc into a nucleotide string */
std::string nucToString(const Nuc);
/** Mutate a nucleotide by an offset. */
char shiftNucleotide(const char, const int);

} /* namespace seqio */

#endif /* SEQIO_H */
