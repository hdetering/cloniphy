#pragma once

#include "../stringio.hpp"
#include "../model/DataFrame.hpp"

#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/path.hpp>
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
  int threads;
  model::DataFrame<double> df_sampling;
  std::vector<SampleConfig> m_vec_samples;
  std::string m_mut_gl_model;
  double m_mut_gl_model_params_kappa;
  std::vector<double> m_mut_gl_model_params_nucfreq;
  
  ConfigStore();
  model::DataFrame<double> getSamplingScheme() const;
  bool parseArgs(int ac, char* av[]);
  template<typename T>
    T parse(const YAML::Node& node);
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

/*--------------------------------*
 * function templates definitions *
 *--------------------------------*/

template<typename T>
T ConfigStore::parse(const YAML::Node& node) {
  T val = node.as<T>();
  return val;
}
/** This specialization can parse numbers in scientific format. */
template<> inline
double ConfigStore::parse<double>(const YAML::Node& node) {
  std::string s = node.as<std::string>();
  double val = stringio::strToDub(s);
  return val;
}

template<typename T>
T ConfigStore::getValue(const char* key) {
  std::vector<std::string> keys = stringio::split(std::string(key), ':');
  YAML::Node node;
  // try to get top-level node
  try {
    node = _config[keys[0]];
  } catch (const std::exception& e) {
    return node.as<T>();
  }
  for (unsigned i=1; i<keys.size(); i++)
    try {
      node = node[keys[i]];
    } catch (const YAML::TypedBadConversion<T>& e) {
      fprintf(stderr, "[WARN] ConfigStore: unknown parameter: '%s'\n", key);
      YAML::Node empty_node;
      return empty_node.as<T>();
    }
/*switch (node.Type()) {
  case YAML::NodeType::Null:
    fprintf(stderr, "'%s': Null\n", keys[keys.size()-1].c_str());
    break;
  case YAML::NodeType::Scalar:
    fprintf(stderr, "'%s': Scalar\n", keys[keys.size()-1].c_str());
    break;
  case YAML::NodeType::Sequence:
    fprintf(stderr, "'%s': Sequence\n", keys[keys.size()-1].c_str());
    break;
  case YAML::NodeType::Map:
    fprintf(stderr, "'%s': Map\n", keys[keys.size()-1].c_str());
    break;
  case YAML::NodeType::Undefined:
    fprintf(stderr, "'%s': Undefined\n", keys[keys.size()-1].c_str());
    break;
  default:
    fprintf(stderr, "'%s': WTF?!\n", keys[keys.size()-1].c_str());
}*/
  return node.as<T>();
}

template<typename T>
T 
ConfigStore::getValue(const std::string key) {
  return getValue<T>(key.c_str());
}

template<typename T>
std::vector<std::vector<T>> 
ConfigStore::getMatrix(const char* key) {
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
std::map<std::string, T> 
ConfigStore::getMap(const char* key) {
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
