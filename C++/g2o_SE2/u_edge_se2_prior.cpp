#include "u_edge_se2_prior.h"

namespace g2o{
    namespace SE2{
        UEdgeSE2Prior::UEdgeSE2Prior() : 
            BaseUnaryEdge< dof_x, Pose, VertexSE2>(){

        }

        void UEdgeSE2Prior::computeError()
        {
            // Get vertex
            const VertexSE2* v1 = static_cast<const VertexSE2*>(_vertices[0]);
            
            // Compute the left-invariant error in the group
            Pose delta = v1->estimate().inverse() * _measurement;
                            
            // Get the coordinates (the twist) of the error (i.e., coordinates of the error in the Lie algebra)
            _error = delta.log().coeffs();
        }

        bool UEdgeSE2Prior::read(std::istream& is)
        {
            // Vector3 p;
            // is >> p[0] >> p[1] ;
            // _measurement = p;
            // for (int i = 0; i < 2; ++i)
            //     for (int j = i; j < 2; ++j) {
            //     is >> information()(i, j);
            //     if (i != j)
            //         information()(j, i) = information()(i, j);
            //     }
            return false;
        }

        bool UEdgeSE2Prior::write(std::ostream& os) const
        {
            // Vector2 p = measurement();
            // os << p(0) << " " << p(1);;
            // for (int i = 0; i < 2; ++i)
            //     for (int j = i; j < 2; ++j)
            //     os << " " << information()(i, j);
            // return os.good();
            return false;
        }
    }
}