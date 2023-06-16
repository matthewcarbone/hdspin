add_library(
    state
    STATIC state.cpp
)

set(PRECISON $ENV{HDSPIN_PRECISON})

if (${PRECISON} MATCHES "^[0-9]+$")
    target_compile_definitions(
        state
        PUBLIC -DPRECISON=${PRECISON}
    )
endif()

set_target_properties(
    state
    PROPERTIES CXX_STANDARD 17
)
