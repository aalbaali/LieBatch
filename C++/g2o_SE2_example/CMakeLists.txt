cmake_minimum_required( VERSION 3.10)
project(example_SE2)

# Main executable
add_executable(example_SE2 main.cpp)

target_link_libraries(example_SE2
  PUBLIC    
    InEKF::SE2
    se2_lib
  )

# Set required C++ standard flag
set_property(TARGET example_SE2 PROPERTY CXX_STANDARD 17)
