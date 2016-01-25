#ifndef SEQIO_H
#define SEQIO_H

#include "stringio.hpp"
#include <iostream>
#include <string>
#include <vector>

/** Handles sequence files (read, write, index). */
namespace seqio {

enum Nuc { A, C, G, T, N };

struct SeqRecord
{
  std::string id;          /** Identifier */
  std::string description; /** Sequence description (everything after first space in ID line) */
  std::string seq;         /** Actual sequence */
  SeqRecord(const std::string&, const std::string&, const std::string&);
};

/** Stores a set of SeqRecords along with indexing information. */
struct Genome
{
  unsigned num_records;
  unsigned length;                        /** total length of all sequences */
  unsigned masked_length;                 /** length of unmasked (non-'N') positions */
  std::vector<SeqRecord> records;
  std::vector<unsigned> vec_start_chr;    /** cumulative start positions of sequences */
  std::vector<unsigned> vec_start_masked; /** cumulative start positions of unmasked regions */
  std::vector<unsigned> vec_cumlen_masked; /** cumulative lengths of unmasked regions */

  Genome(const char*);
  void indexRecords();
};

/** Reads sequences from file. */
std::vector<SeqRecord> readFasta(const char*);
/** Reads sequences from istream. */
std::vector<SeqRecord> readFasta(std::istream&);
/** Writes sequences to ostream. */
int writeFasta(const std::vector<SeqRecord>&, std::ostream&);
/** Generate an index for a FASTA file containing multiple sequences */
void indexFasta(const char*);
/** Identify start positions of chromosomes and masked regions. */
void indexGenome();
/** Convert a nucleotide char into Nuc */
Nuc charToNuc(const char);
/** Convert Nuc into a nucleotide char */
char nucToChar(const Nuc);
/** Mutate a nucleotide by an offset. */
char shiftNucleotide(const char, const int);

} /* namespace seqio */

#endif /* SEQIO_H */
