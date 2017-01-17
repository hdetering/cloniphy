#include <cassert>
#include "stringio.hpp"

using namespace std;

namespace stringio {

vector<string> &split(const string &s, char delim, vector<string> &elems) {
    stringstream ss(s);
    string item;
    while (std::getline(ss, item, delim)) {
        elems.push_back(item);
    }
    return elems;
}

vector<string> split(const string &s, char delim) {
    vector<string> elems;
    split(s, delim, elems);
    return elems;
}

/** Deals with different line endings used on Linux, Mac, Windows platforms */
istream& safeGetline(istream& is, string& t)
{
    t.clear();

    // The characters in the stream are read one-by-one using a std::streambuf.
    // That is faster than reading them one-by-one using the std::istream.
    // Code that uses streambuf this way must be guarded by a sentry object.
    // The sentry object performs various tasks,
    // such as thread synchronization and updating the stream state.

    istream::sentry se(is, true);
    streambuf* sb = is.rdbuf();

    for(;;) {
        int c = sb->sbumpc();
        switch (c) {
        case '\n':
            return is;
        case '\r':
            if(sb->sgetc() == '\n')
                sb->sbumpc();
            return is;
        case EOF:
            // Also handle the case when the last line has no line ending
            if(t.empty())
                is.setstate(ios::eofbit);
            return is;
        default:
            t += (char)c;
        }
    }
}

/** Read a CSV file into a vector of vectors */
unsigned readCSV(vector<vector<string>> &mtx_data, const string filename, char sep) {
  unsigned n_row, n_col = 0;
  ifstream filestream;
  stringio::CSVRow row(sep);

  filestream.open(filename.c_str());
  assert(filestream.good());
  // determine number of columns from first line
  filestream >> row;
  n_col = row.size();

  // read CSV file line by line, splitting cells by separating character
  do {
    unsigned cols = row.size();
    vector<string> vec_row(cols);
    if (cols == 0) { // skip empty lines
      continue;
    }
    else if (cols != n_col) { // issue a warning if #columns differ
      fprintf(stderr, "[WARN] (readCsv) Row #%d has different number of columns (%u) than first row (%u).\n", n_row+1, cols, n_col);
    }
    for (unsigned i=0; i<cols; i++) {
      vec_row[i] = row[i];
    }
    mtx_data.push_back(vec_row);
    n_row++;
  } while (filestream >> row);

  filestream.close();
  return n_row;
}

/** Enable access to CSVRow fields by index. */
string const& CSVRow::operator[](size_t index) const {
    return m_data[index];
}
/** Return number of fields in CSV row. */
size_t CSVRow::size() const {
    return m_data.size();
}
/** Read line from input stream. */
void CSVRow::readNextRow(istream& str) {
    string line;
    getline(str, line);

    stringstream lineStream(line);
    string cell;

    m_data.clear();
    while(getline(lineStream, cell, m_sep))
    {
        m_data.push_back(cell);
    }
}

/** Make data vector accessible */
vector<string> CSVRow::getData() {
  return this->m_data;
}

/** Enable streaming into CSVRow from input stream. */
istream& operator>>(std::istream& str, CSVRow& data) {
    data.readNextRow(str);
    return str;
}

} /* namespace stringio */
