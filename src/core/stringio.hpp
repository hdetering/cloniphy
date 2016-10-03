#ifndef STRINGIO_H
#define STRINGIO_H

#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

namespace stringio {

/** Splits a string by a delimitor into an existing vector */
std::vector<std::string> &split(const std::string&, char, std::vector<std::string>&);
/** Splits a string by a delimiter into an existing vector */
std::vector<std::string> split(const std::string&, char);
/** Reads a line from a stream, dealing with different styles of line endings. */
std::istream& safeGetline(std::istream& is, std::string& t);
/** takes a map and formats it as a table with keys in first column */
template<typename T>
std::string printMatrix(std::map<std::string, std::vector<T>> mtx);

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

} /* namespace stringio */

#endif /*STRINGIO_H */
