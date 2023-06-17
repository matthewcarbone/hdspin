cmake_minimum_required(VERSION 3.6)
project(hdspin)


option(PRECISON "Set the Arbitrary Precision" 256)
option(BUILD_TESTS "Build tests or not" OFF)

# Essentially -Iinc
include_directories(inc)

if (${BUILD_TESTS})
    
    # Find other CMakeLists.txt files in lib and execute
    add_subdirectory(lib/Catch2_v3)

    # Compile the tests, link against required libraries
    add_executable(
        tests
        tests/tests.cpp
        src/energy_mapping.cpp
        src/utils.cpp
    )

    target_compile_definitions(
        tests
        PUBLIC -DPRECISON=${PRECISON}
    )

    set_target_properties(
        tests
        PROPERTIES CXX_STANDARD 17
    )

    target_link_libraries(
        tests
        PRIVATE Catch2::Catch2WithMain
    )

endif()
