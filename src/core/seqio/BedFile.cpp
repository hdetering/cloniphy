#include "BedFile.hpp"

using namespace std;

namespace seqio {

BedFile::BedFile() : 
  m_num_recs(0),
  m_num_feat(0)
{}

BedFile::BedFile(const string& filename) :
  m_num_recs(0),
  m_num_feat(0)
{
  ifstream inputFile;
  inputFile.open(filename, ios::in);
  this->parseBed(inputFile);
  inputFile.close();
}
BedFile::BedFile(const string& filename, const vector<string>& colnames) :
  m_num_recs(0),
  m_num_feat(0)
{
  this->m_feat_names = colnames;

  ifstream inputFile;
  inputFile.open(filename, ios::in);
  this->parseBed(inputFile);
  inputFile.close();
} 

bool
BedFile::parseBed(ifstream& input) {
  string line;

  try {
    stringio::safeGetline(input, line);
    while (input.good()) {
      // skip comment lines
      if (line[0] == '#') continue;
      // split line into fields
      vector<string> row = stringio::split(line, '\t');
      if (row.size()<3) {
        cerr << "[ERROR] Malformed line in BED file (>=3 columns expected)." << endl;
      }
      this->m_vec_seqid.push_back(row[0]);
      this->m_vec_start.push_back(stoul(row[1]));
      this->m_vec_end.push_back(stoul(row[2]));
      // store additional feature columns as string values
      for (size_t i=3; i<row.size(); i++) {
        // add a new feature vector if necessary
        if (this->m_feat_val.size() < (i-2)) {
          this->m_feat_val.push_back(vector<string>());
          this->m_num_feat++;
        }
        this->m_feat_val[i-3].push_back(row[i]);
      }
      this->m_num_recs++;
      stringio::safeGetline(input, line);
    }
  } catch(exception e) {
    cerr << e.what() << endl;
    cerr << "Possibly premature end of BED file?" << endl;
    return false;
  }

  return true;
}

void BedFile::getLoci ( vector<Locus>& out_loci ) {
  out_loci.clear();
  for (size_t i=0; i<this->m_num_recs; i++) {
    out_loci.push_back(Locus(this->m_vec_seqid[i], this->m_vec_start[i], this->m_vec_end[i]));
  }
}

} /* namespace seqio */