#include "types_slam_se2.h"

#include "g2o/core/factory.h"
#include "g2o/stuff/macros.h"

#include <iostream>

namespace g2o {
  namespace SE2 {

  G2O_REGISTER_TYPE_GROUP(custom_slam_2d);

  G2O_REGISTER_TYPE(VERTEX_SE2, VertexSE2);

  G2O_REGISTER_TYPE(BINARY_EDGE_SE2_SE2, BEdgeSE2SE2);

  G2O_REGISTER_TYPE(UNARY_EDGE_SE2_GPS,  UEdgeSE2Gps);

  G2O_REGISTER_TYPE(UNARY_EDGE_SE2_PRIOR,  UEdgeSE2Prior);
  }
} // end namespace
