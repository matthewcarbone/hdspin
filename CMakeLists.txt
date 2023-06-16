cmake_minimum_required(VERSION 3.6)
project(hdspin)

# Essentially -Iinc
include_directories(inc)

# Find other CMakeLists.txt files in lib and execute
add_subdirectory(lib)

# Compile the tests, link against required libraries
add_executable(tests src/tests.cpp)
target_link_libraries(
    tests
    PRIVATE Catch2::Catch2WithMain state_manipulation
)