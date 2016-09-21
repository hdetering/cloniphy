add_subdirectory (core)
add_subdirectory (test)

link_directories(${CMAKE_SOURCE_DIR}/lib)

add_executable (demo main.cpp)
target_link_libraries (demo cloniphy)

# enforce c++11 standard
set_property(TARGET demo PROPERTY CXX_STANDARD 11)
set_property(TARGET demo PROPERTY CXX_STANDARD_REQUIRED ON)

# ----------------------------------------------------------------------------
# Dependencies
# ----------------------------------------------------------------------------

# Search for Threads (prerequisite for ZLIB?)
find_package (Threads)
target_link_libraries (demo ${CMAKE_THREAD_LIBS_INIT})

# Search for zlib as a dependency for SeqAn.
find_package (ZLIB)

# Load the SeqAn module and fail if not found.
#find_package (SeqAn REQUIRED)

# add Boost support
set(Boost_USE_STATIC_LIBS OFF)
set(Boost_USE_MULTITHREADED OFF)
set(Boost_USE_STATIC_RUNTIME OFF)
find_package(Boost 1.58.0 REQUIRED COMPONENTS program_options)

# ----------------------------------------------------------------------------
# Build Setup
# ----------------------------------------------------------------------------

if(Boost_FOUND)
    include_directories(${Boost_INCLUDE_DIRS})
    target_link_libraries(demo ${Boost_LIBRARIES})
endif()

# add zlib support
if (ZLIB_FOUND)
    include_directories(${ZLIB_INCLUDE_DIRS})
    add_definitions(-DSEQAN_HAS_ZLIB)
    target_link_libraries(demo ${ZLIB_LIBRARIES})
endif()
