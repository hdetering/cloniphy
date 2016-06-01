#ifndef SEQIO_H
#define SEQIO_H

#include "random.hpp"
#include "stringio.hpp"
#include <boost/function.hpp>
#include <iostream>
#include <map>
#include <string>
#include <vector>

/** Handles sequence files (read, write, index). */
namespace seqio {

enum Nuc { A, C, G, T, N };

inline short nuc2idx (char nuc) {
  switch (nuc) {
    case 'A': case 'a': return 0;
    case 'C': case 'c': return 1;
    case 'G': case 'g': return 2;
    case 'T': case 't': return 3;
  }
  return -1; // default for unknown chars
}

inline char idx2nuc (short idx) {
  switch (idx) {
    case 0: return 'A';
    case 1: return 'C';
    case 2: return 'G';
    case 3: return 'T';
  }
  return '?'; // default for unknown chars
}

struct SeqRecord
{
  std::string id;          /** identifier */
  std::string description; /** sequence description (everything after first space in ID line) */
  std::string seq;         /** actual sequence */
  std::string id_ref;      /** identifier in reference genome (ploidy) */
  short copy;              /** chromosome copy (1 for haploid) */
  SeqRecord(const std::string, const std::string, const std::string&);
};

/** Represents a genomic location */
struct Locus
{
  unsigned idx_record; // index of genomic sequence
  unsigned start;      // local start position in sequence
  unsigned length;     // length of locus (e.g. =1 for single nucleotide)
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
  std::map<char, std::vector<long> > map_nuc_pos; /** absolute bp positions indexed by nucleotide */
  std::vector<std::vector<long> > nuc_pos;

  Genome();
  Genome(const char*);
  void generate(const unsigned long, const std::vector<double>, RandomNumberGenerator<>&);
  void generate(const unsigned long, const unsigned short, const std::vector<double>, RandomNumberGenerator<>&);
  void indexRecords();
  void duplicate(); // increase ploidy
  /** Get chromosome and local position for global position */
  Locus getLocusByGlobalPos(long) const;
  /** Get absolute coordinates for a relative position in unmasked part of the genome */
  Locus getAbsoluteLocusMasked(double) const;
};

/** Reads sequences from file. */
void readFasta(const char*, std::vector<SeqRecord>&);
/** Reads sequences from istream. */
void readFasta(std::istream&, std::vector<SeqRecord>&);
/** Writes sequences to ostream. */
int writeFasta(const std::vector<SeqRecord>&, std::ostream&, int = 60);
/** Generate an index for a FASTA file containing multiple sequences */
void indexFasta(const char*);
/** Generate random DNA sequence */
unsigned long generateRandomDnaSeq(
  std::string &seq,
  const unsigned long total_len,
  const std::vector<double> nuc_freqs,
  RandomNumberGenerator<> &rng);
/** Simulate allelic dropout events, masking parts of genome as 'N's. */
void simulateADO_old(const std::string, const float, const int, boost::function<double()>&);
/** Simulate allelic dropout events, masking parts of genome as 'N's. */
void simulateADO(const std::string, const unsigned, const float, const int, boost::function<double()>&);
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
