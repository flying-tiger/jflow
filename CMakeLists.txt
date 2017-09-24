cmake_minimum_required(VERSION 3.9.0)

set(CMAKE_DISABLE_SOURCE_CHANGES ON)
set(CMAKE_DISABLE_IN_SOURCE_BUILD ON)

if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -Wall")
endif()
if(CMAKE_CXX_COMPILER_Id STREQUAL "MSVC")
    set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} /W4")
endif()

project(jflow)

#-------------------------------------------------------------------------
# Test Suite
#-------------------------------------------------------------------------
add_executable(test_jflow
    test/test.cpp
    test/test_euler.cpp
    test/test_jflow.cpp
    test/test_structured_grid.cpp
    source/euler.cpp
    source/structured_grid.cpp
)
target_include_directories(test_jflow
  PRIVATE
    source
    extern/catch/single_include
    extern/eigen
)
set_target_properties(test_jflow
  PROPERTIES
    CXX_EXTENSIONS OFF
    CXX_STANDARD_REQUIRED ON
    CXX_STANDARD 17
)


