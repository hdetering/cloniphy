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

class SampleConfig
{
public:
  std::string m_label;
  std::vector<double> m_vec_prevalence;
  SampleConfig(YAML::Node);
};

class ConfigStore
{
public:
  std::vector<SampleConfig> m_vec_samples;
  ConfigStore();
  bool parseArgs(int ac, char* av[]);
  template<typename T>
    T getValue(const char* key);
  template<typename T>
    T getValue(const std::string key);
  template<typename T>
    std::map<std::string, std::vector<T>> getMatrix(const char* key);
  template<typename T>
    std::map<std::string, T> getMap(const char* key);

private:
  YAML::Node _config;
}; /* class ConfigStore */

bool fileExists(std::string filename);

/*---------------------------------*
 * function templates declarations *
 *---------------------------------*/

template<typename T>
T ConfigStore::getValue(const char* key) {
  std::vector<std::string> keys = stringio::split(std::string(key), ':');
  YAML::Node node = _config[keys[0]];
  for (auto i=1; i<keys.size(); i++)
    node = node[keys[i]];
  return node.as<T>();
}

template<typename T>
T ConfigStore::getValue(const std::string key) {
  return getValue<T>(key.c_str());
}

template<typename T>
std::map<std::string, std::vector<T>> ConfigStore::getMatrix(const char* key) {
  std::map<std::string, std::vector<T>> mtx;
  if (!_config[key]) {
    fprintf(stderr, "[ERROR] (ConfigStore::getMatrix) No element found with name '%s'.\n", key);
    return mtx;
  }
  for (auto row : _config[key]) {
    if (row.size() < 2) {
      fprintf(stderr, "[ERROR] (ConfigStore::getMatrix) Too few elements while reading parameter '%s' (key: '%s').\n", row[0].as<std::string>().c_str(), key);
      continue;
    }
    std::vector<T> v = std::vector<T>(row.size()-1);
    for (std::size_t i=1; i<row.size(); i++)
      v[i-1] = row[i].as<T>();
    mtx[row[0].as<std::string>()] = v;
  }
  return mtx;
}

template<typename T>
std::map<std::string, T> ConfigStore::getMap(const char* key) {
  std::map<std::string, T> m;
  if (!_config[key]) {
    fprintf(stderr, "[ERROR] (ConfigStore::getMap) No element found with name '%s'.\n", key);
    return m;
  }
  for (auto row : _config[key]) {
    if (row.size() < 2) {
      fprintf(stderr, "[ERROR] (ConfigStore::getMap) Too few elements while reading parameter '%s' (key: '%s').\n", row[0].as<std::string>().c_str(), key);
      continue;
    }
    T v = row[1].as<T>();
    m[row[0].as<std::string>()] = v;
  }
  return m;
}

} /* namespace config */

#endif /* CONFIGSTORE_H */
