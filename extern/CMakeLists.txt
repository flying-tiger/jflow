# Catch library setup (header only)
set(CATCH_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/catch/single_include)
add_library(catch INTERFACE)
target_include_directories(catch INTERFACE ${CATCH_INCLUDE_DIR})

# Eigen library setup (header only)
set(EIGEN_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/eigen)
add_library(eigen INTERFACE)
target_include_directories(eigen INTERFACE ${EIGEN_INCLUDE_DIR})