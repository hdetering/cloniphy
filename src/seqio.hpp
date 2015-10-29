#ifndef SEQIO_H
#define SEQIO_H

#include <iostream>
#include <string>
#include <vector>

enum Nuc { A, C, G, T, N };

struct SeqRecord
{
  std::string id;  /** Identifier */
  std::string seq; /** Actual sequence */
  //SeqRecord(std::string, std::string);
};

/** Mutations specifiy a modification of a sequence.
 *  TODO: Do I really belong here?
 */
struct Mutation
{
  long absPos;  /** absolute bp position */
  int offset;   /** shifts the ancestral genotype to a new one */

  bool operator< (const Mutation&) const; /** make mutations sortable */
  static std::vector<Mutation> sortByPosition(const std::vector<Mutation>&);
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
};

#endif /* SEQIO_H */
