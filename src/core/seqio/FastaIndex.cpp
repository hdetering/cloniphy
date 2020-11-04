#include "FastaIndex.hpp"

using namespace std;

namespace seqio {

FastaIndex::FastaIndex() : m_num_seqs(0) {}

FastaIndex::FastaIndex(const string& filename) {
  ifstream inputFile;
  inputFile.open(filename, ios::in);
  this->parseFai(inputFile);
  inputFile.close();
} 

bool
FastaIndex::parseFai(ifstream& input) {
  string line;

  try {
    while (input.good()) {
      stringio::safeGetline(input, line);
      if (line.length() == 0) continue;
      vector<string> row = stringio::split(line, '\t');
      if (row.size()<5) {
        cerr << "[ERROR] Malformed line in FASTA index file (>=5 columns expected)." << endl;
      }
      this->m_vec_name.push_back(row[0]);
      this->m_vec_len.push_back(stoul(row[1]));
      this->m_vec_offset.push_back(stoul(row[2]));
      this->m_vec_linebases.push_back(stoul(row[3]));
      this->m_vec_linewidth.push_back(stoul(row[4]));
      this->m_num_seqs++;
    }
  } catch(exception e) {
    cerr << e.what() << endl;
    cerr << "Possibly premature end of FASTA index file?" << endl;
    return false;
  }

  return true;
}

} /* namespace seqio */