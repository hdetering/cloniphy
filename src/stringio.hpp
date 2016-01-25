#ifndef STRINGIO_H
#define STRINGIO_H

#include <fstream>
#include <iostream>
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

} /* namespace stringio */

#endif /*STRINGIO_H */
