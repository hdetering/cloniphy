#include "seqio.hpp"
#include <algorithm>
#include <fstream>
#include <iostream>
#include <stdlib.h> // for system() calls
#include <string>
#include <vector>

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

  try {
    getline(input, header);
    while (input.good()) {
      while (getline(input, line) && line[0]!='>') {
        seq += line;
      }
      std::size_t space_pos = header.find(' ');
      std::string seq_id = header.substr(1, space_pos-1);
      std::string seq_desc = header.substr(space_pos+1);
      SeqRecord rec = {seq_id, seq_desc, seq};
      records.push_back(rec);
      header = line;
      seq = "";
    }
  } catch(std::exception e) {
    std::cerr << e.what() << std::endl;
    std::cerr << "Possibly premature end of FASTA file?" << std::endl;
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

char SeqIO::shiftNucleotide(const char base, const int offset) {
  Nuc nuc_old = charToNuc(base);
  Nuc nuc_new = static_cast<Nuc>((static_cast<int>(nuc_old) + offset) % 4); // TODO: this could be more generic (get rid of the hard-coded 4)
  return nucToChar(nuc_new);
}
