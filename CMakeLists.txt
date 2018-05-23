#include(CMakeToolsHelpers)
cmake_minimum_required(VERSION 3.5.1)
set(CMAKE_LEGACY_CYGWIN_WIN32 0)

# make sure local libraries are added to linker search path
if(COMMAND cmake_policy)
  cmake_policy(SET CMP0003 NEW)
endif(COMMAND cmake_policy)

project(cloniphy CXX)

# pick up build elements from 'src' folder
add_subdirectory(src)

# add yaml-cpp
#add_subdirectory(include/yaml-cpp)

# provides "add_test()" command
enable_testing()

# add_test(run_test)
