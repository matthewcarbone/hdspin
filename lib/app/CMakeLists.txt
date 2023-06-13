add_executable(app main.cpp)
target_link_libraries(app state_manipulation)
set_target_properties(
    app
    PROPERTIES CXX_STANDARD 17
)
