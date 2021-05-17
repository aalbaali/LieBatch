// GPS unary edge
#ifndef B_EDGE_SE2R2_H
#define B_EDGE_SE2R2_H

#include "vertex_se2.h"
// #include "g2o_tutorial_slam_R2_api.h"
#include "g2o/core/base_unary_edge.h"

namespace g2o {
  namespace SE2 {

    /**
     * \brief SE(2) unary edge (GPS measurement)
     */
    // BaseUnaryEdge< dimension = 2, measurementType = Vector2, vertex_1_type = VertexSE2>
    // Measurement: Vector2 : [x; y]
    class UEdgeSE2Gps : public BaseUnaryEdge< 2, Vector2, VertexSE2>
    {
      public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
        UEdgeSE2Gps();

        void computeError();
        
        void setMeasurement(const Vector2& m){
            // Set GPS measurement m = [ x_1; x_2].
            _measurement = m;
        }

        virtual bool read(std::istream& is);
        virtual bool write(std::ostream& os) const;
    };
  }
} // end namespace
#endif
