add_executable(jflow
    main.cpp
)
target_link_libraries(jflow
  PRIVATE
    jflow_lib
    yaml
)
set_target_properties(jflow PROPERTIES ${JFLOW_COMMON_PROPERTIES})
target_compile_options(jflow ${JFLOW_COMMON_OPTIONS})
target_compile_definitions(jflow ${JFLOW_COMMON_DEFINES})