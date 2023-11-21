cmake_minimum_required(VERSION 3.6)
project(hdspin)

option(PRECISON "Set the Arbitrary Precision" 256)
option(BUILD_TESTS "Build tests or not" OFF)
option(SMOKE "Whether or not to use smoke tests" ON)

set (CMAKE_CXX_STANDARD 17)
set (CMAKE_CXX_STANDARD_REQUIRED True)
set (EXTRA_CXX_FLAGS "-O3 -Wall")

include_directories(inc)  # Essentially -Iinc

if (${BUILD_TESTS})

    add_executable(
        tests
        tests/tests.cpp
        src/energy_mapping.cpp
        src/utils.cpp
        src/spin.cpp
        src/obs1.cpp
        src/obs2.cpp
    )

    # Handle the smoke tests
    if (${SMOKE})
        target_compile_definitions(
            tests
            PUBLIC -DPRECISON=${PRECISON} -DSMOKE=1
        )
    else()
        target_compile_definitions(
            tests
            PUBLIC -DPRECISON=${PRECISON} -DSMOKE=0
        )
    endif()

endif()

find_package(MPI REQUIRED)

if(MPI_FOUND)
    set (EXTRA_INCLUDES ${MPI_CXX_INCLUDE_DIRS})
    set (EXTRA_CXX_FLAGS ${EXTRA_CXX_FLAGS} ${MPI_CXX_COMPILE_FLAGS})
    set (EXTRA_LIBS "${EXTRA_LIBS} ${MPI_CXX_LIBRARIES}")
    set (EXTRA_LIBS "${EXTRA_LIBS} ${MPI_CXX_LINK_FLAGS}")
    message ("Top level CMAKE_CXX_FLAGS: " ${CMAKE_CXX_FLAGS})
    message ("Extra CXX flags: " ${EXTRA_CXX_FLAGS})
    message ("Extra includes: " ${EXTRA_INCLUDES})
    message ("Linked libraries are: " ${EXTRA_LIBS})
else()
    message (SEND_ERROR "hdspin requires MPI")
endif()

include_directories(${MPI_CXX_INCLUDE_DIRS})

add_executable(
    hdspin
    src/main.cpp
    src/energy_mapping.cpp
    src/utils.cpp
    src/main_utils.cpp
    src/spin.cpp
    src/obs1.cpp
    src/obs2.cpp
)

target_compile_definitions(hdspin PUBLIC -DPRECISON=${PRECISON})

target_link_libraries(hdspin ${MPI_CXX_LIBRARIES})
