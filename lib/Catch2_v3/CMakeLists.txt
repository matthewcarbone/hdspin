# Build the Catch2 library. This includes compiling a test to make sure
# Catch2_v3 is actually working properly.

project(catch2_testing)
cmake_minimum_required(VERSION 3.6)
add_subdirectory(lib/Catch2)
# add_executable(test_catch2_v3 test_catch2_v3.cpp)
# target_link_libraries(test_catch2_v3 PRIVATE Catch2::Catch2WithMain)
