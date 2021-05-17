#include "vertex_se2.h"

namespace g2o {
  namespace SE2 {

    VertexSE2::VertexSE2() :
      BaseVertex<3, Pose>()
    {        
        // _estimate = Pose();
    }

    bool VertexSE2::read(std::istream& is)
    {
      double x, y, theta;
      is >> x >> y >> theta;
    //   _estimate.fromVector(p);
      _estimate = Pose( x, y, theta);
      return true;
    }

    bool VertexSE2::write(std::ostream& os) const
    {
      Pose X = estimate();
      Eigen::Vector3d p = X.log().coeffs();
      os << p[0] << " " << p[1] << " " << p[2];
      return os.good();
    }
  } // end namespace
} // end namespace
