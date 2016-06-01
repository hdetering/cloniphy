add_subdirectory (core)
add_subdirectory (test)

link_directories(${CMAKE_SOURCE_DIR}/lib)

add_executable (demo main.cpp)
target_link_libraries (demo cloniphy)

# enforce c++11 standard
set_property(TARGET demo PROPERTY CXX_STANDARD 11)
set_property(TARGET demo PROPERTY CXX_STANDARD_REQUIRED ON)

# add Boost support
set(Boost_USE_STATIC_LIBS OFF)
set(Boost_USE_MULTITHREADED OFF)
set(Boost_USE_STATIC_RUNTIME OFF)
find_package(Boost 1.58.0 REQUIRED COMPONENTS program_options)

if(Boost_FOUND)
    include_directories(${Boost_INCLUDE_DIRS})
    target_link_libraries(demo ${Boost_LIBRARIES})
endif()
