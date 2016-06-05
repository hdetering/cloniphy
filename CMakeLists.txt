cmake_minimum_required(VERSION 3.5.1)
set(CMAKE_LEGACY_CYGWIN_WIN32 0)

project(cloniphy CXX)
add_subdirectory(src)

# provides "add_test()" command
enable_testing()

# set compiler flags (dependent on platform)
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU" OR
    "${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
    set(warnings "-W -Wall -pedantic")
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
    set(warnings "/W4 /WX /EHsc")
endif()
set(CMAKE_C_FLAGS   "${CMAKE_C_FLAGS} ${warnings}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${warnings}")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g -O0")

# make external header files available
# NO EFFECT?!
#include_directories(include)

# add_test(run_test)
