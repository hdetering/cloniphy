#ifndef STRINGIO_H
#define STRINGIO_H

#include <cassert>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <map>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

namespace stringio {

/** String formatting. */
template<typename ... Args>
std::string format(const std::string& format, Args ... args);
/** Splits a string by a delimiter into an existing vector */
std::vector<std::string> &split(const std::string&, char, std::vector<std::string>&);
/** Splits a string by a delimiter into a new vector */
std::vector<std::string> split(const std::string&, char);
/** Reads a line from a stream, dealing with different styles of line endings. */
std::istream& safeGetline(std::istream& is, std::string& t);
/** takes a map and formats it as a table with keys in first column */
template<typename T>
std::string printMatrix(std::map<std::string, std::vector<T>> mtx);
/** Read a CSV file into a vector of vectors */
unsigned readCSV (
  std::vector<std::vector<std::string>> &data,
  const std::string filename,
  const char sep=','
);
/** Write vector of vectors to file in CSV format. */
template <typename T>
unsigned writeCSV (
  const std::vector<std::vector<T>> &data,
  const std::string filename,
  const char sep=','
);
/** Write map of named vectors to file in CSV format. */
template <typename T>
unsigned writeCSV (
  const std::map<std::string, std::vector<T>> &data,
  const std::string filename,
  const char sep=','
);



/** Reads a line from a CSV file and splits it into cells.  */
class CSVRow
{
    public:
        /** Default c'tor. Sets separator to ','. */
        CSVRow() : m_sep(',') {}
        /** Inititalize CSVRow with given char as separator. */
        CSVRow(char sep) : m_sep(sep) {}

        /** Enable access to CSVRow fields by index. */
        std::string const& operator[](std::size_t index) const;
        /** Return number of fields in CSV row. */
        std::size_t size() const;
        /** Read line from input stream. */
        void readNextRow(std::istream& str);
        /** Make data vector accessible */
        std::vector<std::string> getData();

    private:
        std::vector<std::string> m_data;
        char m_sep;
};

/** Enable streaming into CSVRow from input stream. */
std::istream& operator>>(std::istream& str, CSVRow& data);


/*---------------------------------*
 * function templates declarations *
 *---------------------------------*/

/** Prints map as a table with keys as first column */
template<typename T>
std::string printMatrix(std::map<std::string, std::vector<T>> mtx) {
  std::string output;
  for (auto record : mtx) {
    output += "  " + record.first;
    for (auto v : record.second)
      output += " | " + std::to_string(v);
    output += "\n";
  }
  return output;
}

/** Write vector of vectors to file in CSV format. */
template <typename T>
unsigned writeCSV (
  const std::vector<std::vector<T>> &data,
  const std::string filename,
  const char sep)
{
  std::ofstream filestream;
  filestream.open(filename.c_str());

  for (auto row : data) {
    auto it = row.begin();
    filestream << *it;
    while (++it != row.end()) {
      filestream << sep << *it;
    }
    filestream << std::endl;
  }

  filestream.close();
}

/** Write map of named vectors to file in CSV format. */
template <typename T>
unsigned writeCSV (
  const std::map<std::string, std::vector<T>> &data,
  const std::string filename,
  const char sep)
{
  std::ofstream filestream;
  filestream.open(filename.c_str());

  for (auto row : data) {
    filestream << row.first;
    for (auto cell : row.second) {
      filestream << sep << cell;
    }
    filestream << std::endl;
  }

  filestream.close();
}

/* Templated function definitions. */

template<typename ... Args>
std::string format( const std::string& format, Args ... args )
{
  size_t size = std::snprintf( nullptr, 0, format.c_str(), args ... ) + 1; // Extra space for '\0'
  std::unique_ptr<char[]> buf( new char[ size ] ); 
  std::snprintf( buf.get(), size, format.c_str(), args ... );
  return std::string( buf.get(), buf.get() + size - 1 ); // We don't want the '\0' inside
}

} /* namespace stringio */

#endif /*STRINGIO_H */
