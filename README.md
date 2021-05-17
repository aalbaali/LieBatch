# In this repo
Implementation of batch inference on robotic states using states that live in Lie groups.

For now, the implementations are on *SE(2)* group with GPS correction. The algorithms are implemented on 
- Matlab (full implementation), and
- C++ using [g2o](https://github.com/RainerKuemmerle/g2o) with custom types that use Lie group objects implemented using [Manif](https://github.com/artivis/manif).

# Problem setup 
- Robot moving on a plane.
- Robot can be described using *SE(2)* elements.
- Robot is measuring
    - angular velocity (interoceptive measurement),
    - linear velocity (interoceptive measurement), and
    - GPS measurement (exteroceptive measurement).
The GPS measurement is coming at different frequency that the simulation frequency.    

# Initialization
The batch optimization are initialized with either
- odometry (i.e., dead-reckoning), or
- left-invariant extended Kalman filter (L-InEKF).

A document is included that explains the batch setup.

# Other repos
- [RandomVariable](https://github.com/aalbaali/RandomVariable)
- [cpp_filter](https://github.com/aalbaali/cpp_filter) (includes the L-InEKF implementation)
- [Manif](https://github.com/artivis/manif)
- [g2o](https://github.com/RainerKuemmerle/g2o)
- [Eigen](https://github.com/libigl/eigen)


