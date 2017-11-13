cmake_minimum_required(VERSION 3.6)
project(jflow)

set(CMAKE_DISABLE_SOURCE_CHANGES ON)
set(CMAKE_DISABLE_IN_SOURCE_BUILD ON)

set(JFLOW_COMMON_PROPERTIES
    CXX_EXTENSIONS OFF
    CXX_STANDARD_REQUIRED ON
    CXX_STANDARD 17
)

# This should eventually go in a toolchain file
if(MSVC)
    set(JFLOW_COMMON_OPTIONS
      PUBLIC   /std:c++latest
      PRIVATE  /W3 /WX
    )
    set(JFLOW_COMMON_DEFINES
      PRIVATE _CRT_SECURE_NO_WARNINGS
      PRIVATE _SCL_SECURE_NO_WARNINGS
    )
else()
    set(JFLOW_COMMON_OPTIONS
      PRIVATE -W -Wall -Werror
    )
    set(JFLOW_COMMON_DEFINES)
endif()

add_subdirectory(extern)
add_subdirectory(jflow_lib)
add_subdirectory(jflow_main)
add_subdirectory(jflow_test)



