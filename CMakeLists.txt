cmake_minimum_required(VERSION 3.9.0)

set(CMAKE_DISABLE_SOURCE_CHANGES ON)
set(CMAKE_DISABLE_IN_SOURCE_BUILD ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_STANDARD_REQUIRED ON)
set(CMAKE_CXX_STANDARD 17)

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


