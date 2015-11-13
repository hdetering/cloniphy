#ifndef SEQIO_H
#define SEQIO_H

#include <iostream>
#include <string>
#include <vector>

enum Nuc { A, C, G, T, N };

struct SeqRecord
{
  std::string id;          /** Identifier */
  std::string description; /** Sequence description (everything after first space in ID line) */
  std::string seq;         /** Actual sequence */
  //SeqRecord(std::string, std::string);
};

/** Reads and writes sequence files. */
class SeqIO {
  public:
    /** Reads sequences from file. */
    static std::vector<SeqRecord> readFasta(const char*);
    /** Reads sequences from istream. */
    static std::vector<SeqRecord> readFasta(std::istream&);
    /** Writes sequences to ostream. */
    static int writeFasta(const std::vector<SeqRecord>&, std::ostream&);
    /** Generate an index for a FASTA file containing multiple sequences */
    static void indexFasta(const char*);
    /** Convert a nucleotide char into Nuc */
    static Nuc charToNuc(const char);
    /** Convert Nuc into a nucleotide char */
    static char nucToChar(const Nuc);
    /** Mutate a nucleotide by an offset. */
    static char shiftNucleotide(const char, const int);
};

#endif /* SEQIO_H */
