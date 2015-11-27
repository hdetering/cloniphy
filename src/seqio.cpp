#include "seqio.hpp"
#include <algorithm>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdlib.h> // for system() calls
#include <string>
#include <vector>

using namespace std;

SeqRecord::SeqRecord(const string& id, const string& desc, const string& seq) : id(id), description(desc), seq(seq) {}

vector<SeqRecord> SeqIO::readFasta(const char *filename) {
  ifstream inputFile;
  inputFile.open(filename, ios::in);
  vector<SeqRecord> records;
  records = readFasta(inputFile);
  inputFile.close();

  return records;
}

vector<SeqRecord> SeqIO::readFasta(istream &input) {
  vector<SeqRecord> records;
  string header;
  string seq;
  string line;

  try {
    getline(input, header);
    while (input.good()) {
      while (getline(input, line) && line[0]!='>') {
        seq += line;
      }
      size_t space_pos = header.find(' ');
      string seq_id = header.substr(1, space_pos-1);
      string seq_desc = header.substr(space_pos+1);
      //SeqRecord rec = {seq_id, seq_desc, seq};
      SeqRecord rec(seq_id, seq_desc, seq);
      records.push_back(rec);
      header = line;
      seq = "";
    }
  } catch(exception e) {
    cerr << e.what() << endl;
    cerr << "Possibly premature end of FASTA file?" << endl;
  }

  return records;
}

int SeqIO::writeFasta(const vector<SeqRecord> &seqs, ostream &output) {
  int recCount = 0;
  for (vector<SeqRecord>::const_iterator rec=seqs.begin(); rec!=seqs.end(); ++rec) {
    output << ">" << rec->id << endl << rec->seq << endl;
    recCount++;
  }
  return recCount;
}

void SeqIO::indexFasta(const char *filename) {
  // TODO: implement this function
  string fn(filename);
  string cmd = "samtools faidx " + fn;
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

vector<string> &SeqIO::split(const string &s, char delim, vector<string> &elems) {
    stringstream ss(s);
    string item;
    while (std::getline(ss, item, delim)) {
        elems.push_back(item);
    }
    return elems;
}

vector<string> SeqIO::split(const string &s, char delim) {
    vector<string> elems;
    split(s, delim, elems);
    return elems;
}
