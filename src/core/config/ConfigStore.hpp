#ifndef CONFIGSTORE_H
#define CONFIGSTORE_H

#include "../stringio.hpp"
#include <boost/format.hpp>
using boost::format;
using boost::str;
#include <boost/program_options.hpp>
#include <iostream>
#include <sys/stat.h>
#include "yaml-cpp/yaml.h"

#define PROGRAM_NAME "CloniPhy 0.1"

namespace config {

class ConfigStore
{
public:
  ConfigStore();
  bool parseArgs(int ac, char* av[]);
  template<typename T>
    T getValue(const char* key);
  template<typename T>
    T getValue(const std::string key);
  template<typename T>
    std::map<std::string, std::vector<T>> getMatrix(const char* key);

private:
  YAML::Node _config;
}; /* class ConfigStore */

bool fileExists(std::string filename);

/*---------------------------------*
 * function templates declarations *
 *---------------------------------*/

template<typename T>
T ConfigStore::getValue(const char* key) {
  return _config[key].as<T>();
}

template<typename T>
T ConfigStore::getValue(const std::string key) {
  return getValue<T>(key.c_str());
}

template<typename T>
std::map<std::string, std::vector<T>> ConfigStore::getMatrix(const char* key) {
  std::map<std::string, std::vector<T>> mtx;
  for (auto row : _config[key]) {
    std::vector<T> v = std::vector<T>(row.size()-1);
    for (std::size_t i=1; i<row.size(); i++)
      v[i-1] = row[i].as<T>();
    mtx[row[0].as<std::string>()] = v;
  }
  return mtx;
}

} /* namespace config */

#endif /* CONFIGSTORE_H */
