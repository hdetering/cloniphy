#ifndef CONFIGSTORE_H
#define CONFIGSTORE_H

#include "../stringio.hpp"
#include "../model/DataFrame.hpp"

#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/path.hpp>
#include <boost/format.hpp>
using boost::format;
using boost::str;
#include <boost/program_options.hpp>
#include <cassert>
#include <iostream>
#include <sys/stat.h>
#include <gitversion/version.h>
#include <yaml-cpp/yaml.h>

#define PROGRAM_NAME "CloniPhy"

namespace config {

/** Simple stand-in for a data matrix with named rows. */
template <typename T>
using TDataMatrix = std::vector<std::vector<T>>;

class SampleConfig
{
public:
  std::string m_label;
  std::vector<double> m_vec_prevalence;
  SampleConfig(std::string label, std::vector<double> weights);
};

class ConfigStore
{
public:
  model::DataFrame<double> df_sampling;
  std::vector<SampleConfig> m_vec_samples;
  
  ConfigStore();
  model::DataFrame<double> getSamplingScheme() const;
  bool parseArgs(int ac, char* av[]);
  template<typename T>
    T getValue(const char* key);
  template<typename T>
    T getValue(const std::string key);
  template<typename T>
    TDataMatrix<T> getMatrix(const char* key);
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
  for (unsigned i=1; i<keys.size(); i++)
    node = node[keys[i]];
  return node.as<T>();
}

template<typename T>
T ConfigStore::getValue(const std::string key) {
  return getValue<T>(key.c_str());
}

template<typename T>
std::vector<std::vector<T>> ConfigStore::getMatrix(const char* key) {
  std::vector<std::vector<T>> mtx;
  if (!_config[key]) {
    fprintf(stderr, "[ERROR] (ConfigStore::getMatrix) No element found with name '%s'.\n", key);
    return mtx;
  }
  YAML::Node node = _config[key];
  assert ( node.Type() == YAML::NodeType::Sequence );
  for (std::size_t i=0; i<node.size(); i++) {
    auto row = node[i];
    std::vector<T> v = std::vector<T>(row.size());
    for (std::size_t j=0; j<row.size(); j++)
      v[j] = row[j].as<T>();
    mtx.push_back(v);
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
  YAML::Node node = _config[key];
  assert(node.Type() == YAML::NodeType::Map);
  for (auto kv : node) {
    std::string k = kv.first.as<std::string>();
    T v = kv.second.as<T>();
    m[k] = v;
  }
  return m;
}

} /* namespace config */

#endif /* CONFIGSTORE_H */
