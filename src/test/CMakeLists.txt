# source code files
file(GLOB SOURCES "*.cpp")

link_directories(${CMAKE_SOURCE_DIR}/lib)

add_executable (run_tests ${SOURCES})
target_link_libraries (run_tests cloniphy)

# enforce c++11 standard
set_property(TARGET run_tests PROPERTY CXX_STANDARD 11)
set_property(TARGET run_tests PROPERTY CXX_STANDARD_REQUIRED ON)

# add Boost support
set(Boost_USE_STATIC_LIBS OFF)
set(Boost_USE_MULTITHREADED OFF)
set(Boost_USE_STATIC_RUNTIME OFF)
find_package(Boost 1.58.0 REQUIRED COMPONENTS program_options system timer)

if(Boost_FOUND)
    include_directories(${Boost_INCLUDE_DIRS})
    target_link_libraries(run_tests ${Boost_LIBRARIES})
endif()

add_custom_command(
    TARGET run_tests POST_BUILD
    COMMAND ln -sf ${CMAKE_SOURCE_DIR}/data ../../
)
