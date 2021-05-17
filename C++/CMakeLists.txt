cmake_minimum_required( VERSION 3.10)
project(Cpp_Batch)

# Append path with Config.cmake dirs/files

# Find the RandomVariable directory
find_package(g2o REQUIRED)
find_package(InEKF)
find_package(RandomVariable)
find_package(manif REQUIRED)
find_package(Eigen3)

FIND_PATH(CSPARSE_INCLUDE_DIR NAMES cs.h
  PATHS
  /usr/include/suitesparse
  /usr/include
  /opt/local/include
  /usr/local/include
  /sw/include
  /usr/include/ufsparse
  /opt/local/include/ufsparse
  /usr/local/include/ufsparse
  /sw/include/ufsparse
)

# Include examples
add_subdirectory("g2o_SE2")