cmake_minimum_required(VERSION 3.6)
project(hdspin)

# Essentially -Iinc
include_directories(inc)

# Find other CMakeLists.txt files in lib and execute
add_subdirectory(lib)

# Compile the tests, link against required libraries
add_executable(
    tests
    src/tests.cpp
    src/state.cpp
)

set(PRECISON $ENV{HDSPIN_PRECISON})

if (${PRECISON} MATCHES "^[0-9]+$")
    target_compile_definitions(
        tests
        PUBLIC -DPRECISON=${PRECISON}
    )
endif()

set_target_properties(
    tests
    PROPERTIES CXX_STANDARD 17
)


target_link_libraries(
    tests
    PRIVATE Catch2::Catch2WithMain
)
