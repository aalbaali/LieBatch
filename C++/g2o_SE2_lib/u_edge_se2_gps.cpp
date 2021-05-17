#include "u_edge_se2_gps.h"

namespace g2o{
    namespace SE2{
        UEdgeSE2Gps::UEdgeSE2Gps() : 
            BaseUnaryEdge< 2, Vector2, VertexSE2>(){

        }

        void UEdgeSE2Gps::computeError()
        {
            // Get vertex
            const VertexSE2* v1 = static_cast<const VertexSE2*>(_vertices[0]);

            // // Compute the left-invariant error (check notes)
            // Vector2 _error = v1->estimate().rotation().transpose() * (
            //     _measurement - v1->estimate().translation()
            // );
            // Compute the left-invariant error (check notes)
            // Vector2 _error = v1->estimate().rotation().transpose() *(_measurement - v1->estimate().translation());
            _error = _measurement - v1->estimate().translation();
        }

        bool UEdgeSE2Gps::read(std::istream& is)
        {
            // Vector2 p;
            // is >> p[0] >> p[1];
            // _measurement = p;
            // for (int i = 0; i < 2; ++i)
            //     for (int j = i; j < 2; ++j) {
            //     is >> information()(i, j);
            //     if (i != j)
            //         information()(j, i) = information()(i, j);
            //     }
            return false;
        }

        bool UEdgeSE2Gps::write(std::ostream& os) const
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