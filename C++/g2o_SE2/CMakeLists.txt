cmake_minimum_required( VERSION 3.10)
project(example_SE2)

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

# Set required C++ standard flag
set_property(TARGET main PROPERTY CXX_STANDARD 17)
