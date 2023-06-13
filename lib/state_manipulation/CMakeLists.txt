add_library(
    state_manipulation
    STATIC state_manipulation.cpp
)
target_compile_definitions(
    state_manipulation
    PUBLIC -DPRECISON=256
)
set_target_properties(
    state_manipulation
    PROPERTIES CXX_STANDARD 17
)
