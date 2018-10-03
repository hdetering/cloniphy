# source code files
file(GLOB SOURCES "*.cpp")

link_directories(${CMAKE_SOURCE_DIR}/lib)

add_executable (run_tests ${SOURCES})
target_link_libraries (run_tests cloniphycore)

# enforce c++11 standard
set_property(TARGET run_tests PROPERTY CXX_STANDARD 11)
set_property(TARGET run_tests PROPERTY CXX_STANDARD_REQUIRED ON)

# add OpenMP support
include(FindOpenMP)
if(OPENMP_FOUND)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif(OPENMP_FOUND)

# ----------------------------------------------------------------------------
# Dependencies
# ----------------------------------------------------------------------------

# add Boost support
set(Boost_USE_STATIC_LIBS OFF)
set(Boost_USE_MULTITHREADED OFF)
set(Boost_USE_STATIC_RUNTIME OFF)
find_package(Boost 1.58.0 REQUIRED COMPONENTS program_options filesystem system timer)

if(Boost_FOUND)
    include_directories(${Boost_INCLUDE_DIRS})
    target_link_libraries(run_tests ${Boost_LIBRARIES})
endif()

# Search for Threads (prerequisite for ZLIB)
find_package (Threads)
target_link_libraries (run_tests ${CMAKE_THREAD_LIBS_INIT})

# add zlib supported
find_package(ZLIB)
if (ZLIB_FOUND)
    include_directories(${ZLIB_INCLUDE_DIRS})
    add_definitions(-DSEQAN_HAS_ZLIB)
    target_link_libraries(run_tests ${ZLIB_LIBRARIES})
endif()


add_custom_command(
    TARGET run_tests POST_BUILD
    COMMAND ln -sf ${CMAKE_SOURCE_DIR}/data ../../
)
