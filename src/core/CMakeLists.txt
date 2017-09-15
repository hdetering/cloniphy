# source code files
file(GLOB_RECURSE SOURCES "*.cpp")
# header files
file(GLOB_RECURSE HEADERS "*.hpp")

add_library (cloniphy ${SOURCES} ${HEADERS})

# enforce c++11 standard
set_property(TARGET cloniphy PROPERTY CXX_STANDARD 11)
set_property(TARGET cloniphy PROPERTY CXX_STANDARD_REQUIRED ON)

target_include_directories(cloniphy PUBLIC ${CMAKE_SOURCE_DIR}/include bamio config)

# add yaml-cpp library
target_link_libraries(cloniphy yaml-cpp)
