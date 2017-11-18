#-----------------------------------------------------------------------
# Blaze
#-----------------------------------------------------------------------
add_library(blaze INTERFACE)
target_include_directories(blaze
    INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/blaze
)

#-----------------------------------------------------------------------
# Catch
#-----------------------------------------------------------------------
add_library(catch INTERFACE)
target_include_directories(catch
    INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/catch/single_include
)

#-----------------------------------------------------------------------
# YAML
#-----------------------------------------------------------------------
add_subdirectory(yaml-cpp EXCLUDE_FROM_ALL)
add_library(yaml INTERFACE)
target_include_directories(yaml
    INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/yaml-cpp/include
)
target_link_libraries(yaml
    INTERFACE yaml-cpp
)

#-----------------------------------------------------------------------
# NOTES
#-----------------------------------------------------------------------
# [1] Catch does not expose a library target, so we need to make one.
#     This follows the recommendation of the Catch library authors:
#     https://github.com/philsquared/Catch/blob/master/docs/build-systems.md
#
# [2] Yaml-cpp was annoying. It requires a build step, but is trying to
#     target really old CMAKE versions, so it doesn't define its include
#     directories as part of its targets. So, we define a wrapper target
#     that DOES define the include dir and requires linking to the
#     yaml-cpp target.
#
# [3] Blaze actually does have a pretty modern CMakeLists, but I'm having
#     trouble doing a proper install because it can't find LAPACK on
#     windows. Using as header only for now. Will need to update this
#     eventually....



