#ifndef SEQIO_H
#define SEQIO_H

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

/** Reads sequences from file. */
std::vector<SeqRecord> readFasta(const char*);
/** Reads sequences from istream. */
std::vector<SeqRecord> readFasta(std::istream&);
/** Writes sequences to ostream. */
int writeFasta(const std::vector<SeqRecord>&, std::ostream&);
/** Generate an index for a FASTA file containing multiple sequences */
void indexFasta(const char*);
/** Convert a nucleotide char into Nuc */
Nuc charToNuc(const char);
/** Convert Nuc into a nucleotide char */
char nucToChar(const Nuc);
/** Mutate a nucleotide by an offset. */
char shiftNucleotide(const char, const int);
/** Splits a string by a delimitor into an existing vector */
std::vector<std::string> &split(const std::string&, char, std::vector<std::string>&);
/** Splits a string by a delimiter into an existing vector */
std::vector<std::string> split(const std::string&, char);

}

#endif /* SEQIO_H */
