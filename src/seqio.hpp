#ifndef SEQIO_H
#define SEQIO_H

#include "stringio.hpp"
#include <boost/function.hpp>
#include <iostream>
#include <string>
#include <vector>

/** Handles sequence files (read, write, index). */
namespace seqio {

enum Nuc { A, C, G, T, N };

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

  Genome(const char*);
  void indexRecords();
  void duplicate(); // increase ploidy
  /** Get absolute coordinates for a relative position in unmasked pard of the genome */
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
/** Simulate allelic dropout events, masking parts of genome as 'N's. */
void simulateADO_old(const std::string, const float, const int, boost::function<float()>&);
/** Simulate allelic dropout events, masking parts of genome as 'N's. */
void simulateADO(const std::string, const unsigned, const float, const int, boost::function<float()>&);
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
