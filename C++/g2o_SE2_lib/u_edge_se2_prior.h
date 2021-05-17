// GPS unary edge
#ifndef B_EDGE_SE2_PRIOR_H
#define B_EDGE_SE2_PRIOR_H

// Use the custom SE2 manif typedefs defined in `inekf_se2.h`.
#include "inekf_se2.h"
#include "manif/SE2.h"
#include "vertex_se2.h"

// #include "g2o_tutorial_slam_R2_api.h"
#include "g2o/core/base_unary_edge.h"

namespace g2o {
  namespace SE2 {

    /**
     * \brief SE(2) unary edge (SE(2) prior)
     */
    // BaseUnaryEdge< dimension = 3, measurementType = Pose, vertex_1_type = VertexSE2>
    // Measurement: Vector2 : [x; y]
    class UEdgeSE2Prior : public BaseUnaryEdge< dof_x, Pose, VertexSE2>
    {
      public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
        UEdgeSE2Prior();

        void computeError();
        
        void setMeasurement(const Pose& m){
            // Set SE(2) prior
            _measurement = m;
        }

        virtual bool read(std::istream& is);
        virtual bool write(std::ostream& os) const;
    };
  }
} // end namespace
#endif
