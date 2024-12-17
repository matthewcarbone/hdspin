cmake_minimum_required(VERSION 3.6)
project(hdspin)

option(PRECISON "Set the Arbitrary Precision" 256)
option(BUILD_TESTS "Build tests or not" OFF)
option(SMOKE "Whether or not to use smoke tests" ON)
option(VERSION_STAMP "Whether or not to get the current commit hash as version" ON)

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


# Build the core package

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
    src/processing_utils.cpp
    src/spin.cpp
    src/obs1.cpp
    src/obs2.cpp
    src/emax.cpp
)

if (${VERSION_STAMP})
    execute_process(
        COMMAND git rev-parse HEAD
        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
        OUTPUT_VARIABLE GIT_COMMIT_HASH
        OUTPUT_STRIP_TRAILING_WHITESPACE
    )
    target_compile_definitions(hdspin PUBLIC -DPRECISON=${PRECISON} -DGIT_COMMIT_HASH=\"${GIT_COMMIT_HASH}\")
else()
    target_compile_definitions(hdspin PUBLIC -DPRECISON=${PRECISON})

endif()

target_link_libraries(hdspin ${MPI_CXX_LIBRARIES})


# Build the postprocessing executable
# Useful for if the main job doesn't actually finish, but we want to just
# process what we have

add_executable(
    hdspin-postprocess
    src/main_postprocess.cpp
    src/utils.cpp
    src/processing_utils.cpp
)

target_link_libraries(hdspin-postprocess ${MPI_CXX_LIBRARIES})

install(TARGETS hdspin hdspin-postprocess RUNTIME DESTINATION bin)
