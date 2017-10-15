add_library(jflow_lib
    euler.cpp
    jflow.cpp
    structured_grid.cpp
)
target_link_libraries(jflow_lib
  PUBLIC
    eigen
)
target_include_directories(jflow_lib PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
set_target_properties(jflow_lib PROPERTIES ${JFLOW_COMMON_PROPERTIES})
target_compile_options(jflow_lib ${JFLOW_COMMON_OPTIONS})
target_compile_definitions(jflow_lib ${JFLOW_COMMON_DEFINES})