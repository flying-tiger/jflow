cmake_minimum_required(VERSION 3.9.0)

project(jflow)

include_directories(
    source
    extern/catch/single_include
)

add_executable(test_jflow
    test/test.cpp 
    test/test_structured_grid.cpp
    source/structured_grid.cpp
)


