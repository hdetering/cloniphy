#include "seqio.hpp"
#include <algorithm>
#include <fstream>
#include <iostream>
#include <stdlib.h> // for system() calls
#include <string>
#include <vector>

//SeqRecord::SeqRecord(std::)

bool Mutation::operator< (const Mutation &other) const {
  return absPos < other.absPos;
}

std::vector<Mutation> Mutation::sortByPosition(const std::vector<Mutation> &mutations) {
  std::vector<Mutation> mutationsCopy = mutations;
  std::sort(mutationsCopy.begin(), mutationsCopy.end());
  return mutationsCopy;
}

std::vector<SeqRecord> SeqIO::readFasta(const char *filename) {
  std::ifstream inputFile;
  inputFile.open(filename, std::ios::in);
  std::vector<SeqRecord> records;
  records = readFasta(inputFile);
  inputFile.close();

  return records;
}

std::vector<SeqRecord> SeqIO::readFasta(std::istream &input) {
  std::vector<SeqRecord> records;
  std::string header;
  std::string seq;
  std::string line;

  while (input.good()) {
    try {
      getline(input, header);
      while (getline(input, line) && line[0]!='>') {
        seq += line;
      }
    } catch(std::exception e) {
      std::cerr << e.what() << std::endl;
      std::cerr << "Possibly premature end of FASTA file?" << std::endl;
    }
    SeqRecord rec = {header.substr(1), seq};
    records.push_back(rec);
  }

  return records;
}

int SeqIO::writeFasta(const std::vector<SeqRecord> &seqs, std::ostream &output) {
  int recCount = 0;
  for (std::vector<SeqRecord>::const_iterator rec=seqs.begin(); rec!=seqs.end(); ++rec) {
    output << ">" << rec->id << std::endl << rec->seq << std::endl;
    recCount++;
  }
  return recCount;
}

void SeqIO::indexFasta(const char *filename) {
  // TODO: implement this function
  std::string fn(filename);
  std::string cmd = "samtools faidx " + fn;
  system(cmd.c_str());
}


Nuc SeqIO::charToNuc(const char c) {
  switch (c) {
    case 'A':
      return A;
    case 'C':
      return C;
    case 'G':
      return G;
    case 'T':
      return T;
    default:
      return N;
  }
}

char SeqIO::nucToChar(const Nuc n) {
  switch (n) {
    case A:
      return 'A';
    case C:
      return 'C';
    case G:
      return 'G';
    case T:
      return 'T';
    default:
      return 'N';
  }
}
