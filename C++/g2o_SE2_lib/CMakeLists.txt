cmake_minimum_required( VERSION 3.10)
project(library_SE2)

# Add the SE2 library
add_library(se2_lib STATIC
            vertex_se2.cpp vertex_se2.h
            b_edge_se2se2.cpp b_edge_se2se2.h  
            u_edge_se2_gps.cpp  u_edge_se2_gps.h
            u_edge_se2_prior.cpp  u_edge_se2_prior.h
            types_slam_se2.h  types_slam_se2.cpp)

# Link the InEKF::SE2 in order to use the definitions in the header            
target_link_libraries(se2_lib
    PUBLIC
    InEKF::SE2 MANIF::manif
    g2o::core g2o::stuff g2o::solver_eigen g2o::solver_csparse g2o::csparse_extension
)        
target_include_directories(se2_lib 
  PUBLIC
  ${CSPARSE_INCLUDE_DIR}
)