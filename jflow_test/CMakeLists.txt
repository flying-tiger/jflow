add_executable(jflow_test
    test.cpp
    test_common.cpp
    test_euler.cpp
    test_finite_volume.cpp
    test_grid_util.cpp
    test_integrators.cpp
    test_structured_grid.cpp
)
target_link_libraries(jflow_test
  PRIVATE
    blaze
    catch
    jflow_lib
)
set_target_properties(jflow_test PROPERTIES ${JFLOW_COMMON_PROPERTIES})
target_compile_options(jflow_test ${JFLOW_COMMON_OPTIONS})
target_compile_definitions(jflow_test ${JFLOW_COMMON_DEFINES})
