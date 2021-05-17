// SE2 process model edge
#ifndef B_EDGE_SE2SE2_H
#define B_EDGE_SE2SE2_H

#include "vertex_se2.h"
// #include "g2o_tutorial_slam_R2_api.h"
#include "g2o/core/base_binary_edge.h"

// Use the typedefs
#include "inekf_se2.h"

namespace g2o {
  namespace SE2 {

    /**
     * \brief SE(2) edge between two VertexSE2 (i.e., odometry constraint)
     */
    // BaseBinaryEdge< dimension = 3, measurementType = Vector3, vertex_1_type = VertexSE2, vertex_2_type = VertexSE2>
    // Measurement: Vector3 : [x; y; theta]
    class BEdgeSE2SE2 : public BaseBinaryEdge< 3, Vector3, VertexSE2, VertexSE2>
    {
      public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
        BEdgeSE2SE2();

        void computeError();
        
        // Set sampling period dt
        void setDt(const double& dt){ _dt = dt;};
        
        // Get dt
        double dt() const { return _dt;};

        void setMeasurement(const Vector3& m){
            // Set interoceptive measurement m = [ v_1; v2; omega], where v is the linear velocity

            // Ensure that dt is set (it shouldn't be -1)
            assert(_dt != -1);
            _measurement = m;
            _Xi = LieAlg( _dt * _measurement).exp();
        }

        virtual bool read(std::istream& is);
        virtual bool write(std::ostream& os) const;

      protected:
        // Xi = exp( wedge( dt_km1 * u_km1))
        Pose _Xi;

        // Sampling period
        double  _dt = -1; 
    };
  }
} // end namespace
#endif
