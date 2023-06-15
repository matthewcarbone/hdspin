add_library(
    state_manipulation
    STATIC state_manipulation.cpp
)

set(PRECISON $ENV{HDSPIN_PRECISON})

if (${PRECISON} MATCHES "^[0-9]+$")
    target_compile_definitions(
        state_manipulation
        PUBLIC -DPRECISON=${PRECISON}
    )
endif()

set_target_properties(
    state_manipulation
    PROPERTIES CXX_STANDARD 17
)
