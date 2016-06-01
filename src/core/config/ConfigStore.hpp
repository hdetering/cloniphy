#ifndef CONFIGSTORE_H
#define CONFIGSTORE_H

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

private:
  YAML::Node _config;
};

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

} /* namespace config */

#endif /* CONFIGSTORE_H */
