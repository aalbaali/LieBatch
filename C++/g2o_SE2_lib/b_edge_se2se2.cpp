#include "b_edge_se2se2.h"

namespace g2o{
    namespace SE2{
        BEdgeSE2SE2::BEdgeSE2SE2() : 
            BaseBinaryEdge< 3, Vector3, VertexSE2, VertexSE2>(){

        }

        // Implementation of the left-invariant error
        void BEdgeSE2SE2::computeError()
        {
            // Left-invariant error
            const VertexSE2* v1 = static_cast<const VertexSE2*>(_vertices[0]);
            const VertexSE2* v2 = static_cast<const VertexSE2*>(_vertices[1]);

            // Compute the left-invariant error in the group
            Pose delta = v2->estimate().inverse()
                            * v1->estimate()
                            * _Xi;
            
            // Get the coordinates (the twist) of the error (i.e., coordinates of the error in the Lie algebra)
            _error = delta.log().coeffs();
        }

        // Implementation of the Jacobian on the left-invariant error
        //      J_e_xi1 = -SE2.computeJRightInv( _error.exp()) * SE2.adjoint(inv(_Xi)),
        //      J_e_xi2 = SE2.computeJLeftInv( _error.exp()),.
        void BEdgeSE2SE2::linearizeOplus(){
            // Construct the error in Lie algebra object
            LieAlg e_alg = LieAlg( _error);
                        
            // Compute Jacobian w.r.t. the first pose
            _jacobianOplusXi = - e_alg.rjacinv() * _Xi.inverse().adj();            
            _jacobianOplusXj = e_alg.ljacinv();            
        }

        bool BEdgeSE2SE2::read(std::istream& is)
        {
            Vector3 p;
            is >> p[0] >> p[1] >> p[2];
            _measurement = p;
            _Xi = LieAlg( _dt * _measurement).exp();
            for (int i = 0; i < 3; ++i)
                for (int j = i; j < 3; ++j) {
                is >> information()(i, j);
                if (i != j)
                    information()(j, i) = information()(i, j);
                }
            return true;
        }

        bool BEdgeSE2SE2::write(std::ostream& os) const
        {
            Vector3 p = measurement();
            os << p(0) << " " << p(1) << " " << p(2);
            for (int i = 0; i < 3; ++i)
                for (int j = i; j < 3; ++j)
                os << " " << information()(i, j);
            return os.good();
        }
    }
}