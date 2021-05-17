#ifndef VERTEX_SE2_H
#define VERTEX_SE2_H

#include "g2o/core/base_vertex.h"
#include "g2o/core/hyper_graph_action.h"

// Use the custom SE2 manif typedefs defined in `inekf_se2.h`.
#include "inekf_se2.h"
#include "manif/SE2.h"

namespace g2o {
  namespace SE2 {

    // ********************************************
    // SE(2) vertex
    // ********************************************
    class VertexSE2 : public BaseVertex<3, Pose>
    {
      public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
        VertexSE2();

        virtual void setToOriginImpl() {
            _estimate.setIdentity();
        }

        virtual void oplusImpl(const double* update)
        {
          // Lie Algebra element. Note that the update is [x; y; theta].
          LieAlg dx( update[0], update[1], update[2]);

          _estimate += dx;          
        }

        virtual bool read(std::istream& is);
        virtual bool write(std::ostream& os) const;

        void setTime(double time_in){
          _time = time_in;
        }
        const double time() const{
          // Returns a reference to the 'time' variable.
          return _time;
        } 
      private:
        double _time;
    };

  } // end namespace
} // end namespace

#endif
