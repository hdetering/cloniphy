add_subdirectory (core)
add_subdirectory (test)

link_directories(${CMAKE_SOURCE_DIR}/lib)

add_executable (cloniphy main.cpp)
target_link_libraries (cloniphy cloniphycore)

# enforce c++11 standard
set_property(TARGET cloniphy PROPERTY CXX_STANDARD 11)
set_property(TARGET cloniphy PROPERTY CXX_STANDARD_REQUIRED ON)
# NOTE: This should help print string values in GDB sessions, comment out for prod!
#add_definitions(-D_GLIBCXX_USE_CXX11_ABI=0)

# ----------------------------------------------------------------------------
# Dependencies
# ----------------------------------------------------------------------------

# Search for Threads (prerequisite for ZLIB?)
find_package (Threads)
target_link_libraries (cloniphy ${CMAKE_THREAD_LIBS_INIT})

# Search for zlib as a dependency for SeqAn.
find_package (ZLIB)
# add zlib support
if (ZLIB_FOUND)
    include_directories(${ZLIB_INCLUDE_DIRS})
    add_definitions(-DSEQAN_HAS_ZLIB)
    target_link_libraries(cloniphy ${ZLIB_LIBRARIES})
endif()

# add Boost support
set(Boost_USE_STATIC_LIBS OFF)
set(Boost_USE_MULTITHREADED OFF)
set(Boost_USE_STATIC_RUNTIME OFF)
find_package(Boost 1.58.0 REQUIRED COMPONENTS program_options filesystem)
if(Boost_FOUND)
    include_directories(${Boost_INCLUDE_DIRS})
    target_link_libraries(cloniphy ${Boost_LIBRARIES})
endif()

# Load the SeqAn module and fail if not found.
#find_package (SeqAn REQUIRED)

# ----------------------------------------------------------------------------
# Build Setup
# ----------------------------------------------------------------------------

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

# ----------------------------------------------------------------------------
# Manage resource files
# ----------------------------------------------------------------------------

add_custom_command(TARGET cloniphy POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_directory
        ${CMAKE_SOURCE_DIR}/resources
        $<TARGET_FILE_DIR:cloniphy>/resources)