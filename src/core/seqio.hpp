#ifndef SEQIO_H
#define SEQIO_H

#include "random.hpp"
#include "stringio.hpp"
#include <algorithm>
#include <cassert>
#include <functional>
#include <iostream>
#include <map>
#include <string>
#include <vector>

/** Handles sequence files (read, write, index). */
namespace seqio {

enum Nuc { A, C, G, T, N };

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

struct SeqRecord
{
  std::string id;          /** identifier */
  std::string description; /** sequence description (everything after first space in ID line) */
  std::string seq;         /** actual sequence */
  std::string id_ref;      /** identifier in reference genome (ploidy) */
  short chr_copy;          /** chromosome copy (0 for haploid) */
  SeqRecord(const std::string, const std::string, const std::string&);
};

/** Represents a genomic location */
struct Locus
{
  unsigned    idx_record; // index of genomic sequence
  std::string id_ref;     // seq id in ref genome
  unsigned    start;      // local start position in sequence (0-based)
  unsigned    length;     // length of locus (e.g. =1 for single nucleotide)
};

/** Stores a set of SeqRecords along with indexing information. */
struct Genome
{
  unsigned num_records;
  unsigned length;                        /** total length of all sequences */
  unsigned masked_length;                 /** length of unmasked (non-'N') positions */
  // TODO: make ploidy a feature of each sequence (?)
  short ploidy;                           /** number of copies for each chromosome */
  std::vector<SeqRecord> records;

  std::vector<unsigned> vec_start_chr;    /** cumulative start positions of sequences */
  std::vector<unsigned> vec_start_masked; /** cumulative start positions of unmasked regions */
  std::vector<unsigned> vec_cumlen_masked; /** cumulative lengths of unmasked regions */
  double nuc_freq[4];                     /** nucleotide frequencies */
  /** absolute bp positions indexed by nucleotide */
  //std::map<char, std::vector<long> > map_nuc_pos;
  std::vector<std::vector<long> > nuc_pos;
  /** absolute bp positions indexed by tri-nucleotides */
  std::map<std::string, std::vector<long> > map_3mer_pos;

  Genome();
  Genome(const char*);
  /** Simulate DNA seq of given length and nuc freqs. */
  void generate (
    const unsigned long,
    const std::vector<double>,
    RandomNumberGenerator<>&
  );

  /**
   * Generate random reference genome based on nucleotide frequencies.
   *
   * \param total genome length including all sequences
   * \param number of chromosomes/sequences to generate
   * \param nucleotide frequencies (A,C,G,T)
   * \param object to generate random numbers
   */
  void generate (
    const unsigned long,
    const unsigned short,
    const std::vector<double>,
    RandomNumberGenerator<>&
  );

  /**
    * Generate random genome with given number of fragments, mean len, sd len, nuc freqs.
    *
    * \param num_frags  number of chromosomes/sequences to generate
    * \param mean_len   mean sequence length
    * \param sd_len     standard deviation of sequence length
    * \param nuc_freqs  nucleotide frequencies (A,C,G,T)
    * \param rng object to generate random numbers
    */
  void generate (
    const unsigned num_frags,
    const unsigned long mean_len,
    const unsigned long sd_len,
    const std::vector<double> nuc_freqs,
    RandomNumberGenerator<>& rng
  );

  /**
   * Index genome for easier access to sequences and individual nucleotides.
   *
   * Indices are generated to translate between global and local coordinates
   * as well as to quickly locate
   *  - unmasked (non-"N") regions
   *  - positions having a given nucleotide
   *  - positions having a given tri-nucleotide
   */
  void indexRecords();

  /**
   * Increase the ploidy of the whole genome by a factor of two.
   */
  void duplicate();

  /**
   * Get chromosome and local position for global position
   */
  Locus getLocusByGlobalPos(long) const;

  /**
   * Get absolute coordinates for a relative position in unmasked part of the genome
   */
  Locus getAbsoluteLocusMasked(double) const;
};

/** Reads sequences from file. */
void readFasta(const char*, std::vector<SeqRecord>&);
/** Reads sequences from istream. */
void readFasta(std::istream&, std::vector<SeqRecord>&);
/** Writes sequences to file. */
int writeFasta(const std::vector<SeqRecord>&, const std::string fn, int len_line = 60);
/** Writes sequences to ostream. */
int writeFasta(const std::vector<SeqRecord>&, std::ostream& os, int len_line = 60);
/** Generate an index for a FASTA file containing multiple sequences */
void indexFasta(const char*);
/** Generate random DNA sequence */
unsigned long generateRandomDnaSeq(
  std::string &seq,
  const unsigned long total_len,
  const std::vector<double> nuc_freqs,
  RandomNumberGenerator<> &rng);
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
