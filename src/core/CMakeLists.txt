# source code files
file(GLOB_RECURSE SOURCES "*.cpp")
# header files
file(GLOB_RECURSE HEADERS "*.hpp")

add_library (cloniphycore ${SOURCES} ${HEADERS})

# enforce c++11 standard
set_property(TARGET cloniphycore PROPERTY CXX_STANDARD 11)
set_property(TARGET cloniphycore PROPERTY CXX_STANDARD_REQUIRED ON)

# add OpenMP support
include(FindOpenMP)
if(OPENMP_FOUND)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif(OPENMP_FOUND)

target_include_directories(cloniphycore PUBLIC ${CMAKE_SOURCE_DIR}/include bamio config model)

# add yaml-cpp library
target_include_directories(cloniphycore PUBLIC ${CMAKE_SOURCE_DIR}/include/yaml-cpp/include)
target_link_libraries(cloniphycore yaml-cpp)

# use gitversion to make Git version info available in code
# https://github.com/smessmer/gitversion/tree/develop
INCLUDE(${CMAKE_SOURCE_DIR}/include/gitversion/cmake.cmake)
TARGET_GIT_VERSION_INIT(cloniphycore)