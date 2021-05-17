cmake_minimum_required( VERSION 3.10)
project(example_SE2)

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

# Add the SE2 library
add_library(example_se2_lib STATIC
            vertex_se2.cpp vertex_se2.h
            b_edge_se2se2.cpp b_edge_se2se2.h  
            u_edge_se2_gps.cpp  u_edge_se2_gps.h
            u_edge_se2_prior.cpp  u_edge_se2_prior.h
            types_slam_se2.h  types_slam_se2.cpp)

# Link the InEKF::SE2 in order to use the definitions in the header            
target_link_libraries(example_se2_lib
    PUBLIC
    InEKF::SE2 MANIF::manif
    g2o::core g2o::stuff g2o::solver_eigen g2o::solver_csparse g2o::csparse_extension
)        
target_include_directories(example_se2_lib 
  PUBLIC
  ${CSPARSE_INCLUDE_DIR}
)


# Main executable
add_executable(main main.cpp)

target_link_libraries(main
  PUBLIC    
    InEKF::SE2
    example_se2_lib
  )

  # target_include_directories(main PUBLIC 
  # ${CSPARSE_INCLUDE_DIR}
  # )

# Set required C++ standard flag
set_property(TARGET main PROPERTY CXX_STANDARD 17)
