#include <iostream>

#include "inekf_se2.h"

#include "types_slam_se2.h"

#include "g2o/core/sparse_optimizer.h"
#include "g2o/core/sparse_block_matrix.h"
#include "g2o/core/block_solver.h"
#include "g2o/core/factory.h"
#include "g2o/core/optimization_algorithm_factory.h"
#include "g2o/core/optimization_algorithm_gauss_newton.h"
#include <g2o/core/optimization_algorithm_levenberg.h>

#include "g2o/solvers/eigen/linear_solver_eigen.h"
#include <g2o/solvers/csparse/linear_solver_csparse.h>
#include <g2o/solvers/cholmod/linear_solver_cholmod.h>

#include "g2o/core/marginal_covariance_cholesky.h"


int main(int argc, const char* argv[]){
    // Read config.yml. Specify the arguments in the settings.json file. For example,
    //  {
    //     "cmake.debugConfig": {
    //         "args": [
    //             "/home/aa/Documents/Data/Data_generator/SE2/config.yml"
    //         ]
    //     }
    // }
    YAML::Node config;
    std::string filename_config;
    filename_config = argv[1];
#ifndef NDEBUG    
    std::cout << argc << " arguments passed: " << filename_config << std::endl;
#endif    
    config = YAML::LoadFile( filename_config);
    // Set verbosity
    const bool verb = config["verb"].as<bool>();
    // Get output file name
    const std::string filename_kf   = config["filename_kf"].as<std::string>();

    // ********************************************************
    // Load data
#ifndef NDEBUG
    std::cout << "Loading data" << std::endl;
#endif
    // Read sensor files
    // Prior
    const std::string filename_prior = config["filename_prior"].as<std::string>();
    //  Gyro
    const std::string filename_gyro  = config["filename_gyro"].as<std::string>();
    //  Velocity
    const std::string filename_vel   = config["filename_vel"].as<std::string>();
    //  GPS
    const std::string filename_gps   = config["filename_gps"].as<std::string>();
    // Estimated states
    // const std::string filename_out   = config["filename_out"].as<std::string>();

    // Import data
    //  Prior
    PoseEstimate meas_prior = RV::IO::import< PoseEstimate>( filename_prior)[0];
    //  Gyro
    std::vector< MeasGyro> meas_gyro      = RV::IO::import< MeasGyro>( filename_gyro);
    //  Velocity
    std::vector< MeasVel> meas_vel        = RV::IO::import< MeasVel>( filename_vel);
    //  GPS
    std::vector< MeasGps> meas_gps        = RV::IO::import< MeasGps>( filename_gps);
    
    // ********************************************************
    // Run L-InEKF    
    if(verb)
        std::cout << "Running L-InEKF" << std::endl;
    
    // Note that the estimates PoseEstimate is of type RandomVariable
    std::vector< PoseEstimate> X_kf_rv = GetSe2InekfEstimates( 
            meas_prior,
            meas_gyro,
            meas_vel,
            meas_gps
            );

    // Number of poses
    const int K = X_kf_rv.size();

    // dt_func returns sampling period at index k: dt_k = t_k - t_{k-1}
    auto dt_func = [&meas_gyro](int k){
        return meas_gyro[k+1].time() - meas_gyro[k].time();
    };

    // Compute vector of Pose (SE2d) elements
    std::vector< Pose> X_kf( K);
    for( int k = 0; k < K; k++){
        double x, y, real, imag;
        x    = X_kf_rv[k].mean()(0, 2);
        y    = X_kf_rv[k].mean()(1, 2);  
        real = X_kf_rv[k].mean()(0, 0);  
        imag = X_kf_rv[k].mean()(1, 0);  
        X_kf[k] = Pose( x, y, real, imag);
    }

    // ********************************************************
    // creating the optimization problem
    typedef g2o::BlockSolver< g2o::BlockSolverTraits<3, -1> >      SlamBlockSolver;
    typedef g2o::LinearSolverEigen<SlamBlockSolver::PoseMatrixType> SlamLinearSolver;
    typedef g2o::LinearSolverCSparse<SlamBlockSolver::PoseMatrixType> CSparseLinearSolver;
    typedef g2o::LinearSolverCholmod<SlamBlockSolver::PoseMatrixType> CholmodLinearSolver;  // Can be used to compute marginals

    // allocating the optimizer
    g2o::SparseOptimizer optimizer;
    // auto linearSolver = g2o::make_unique<SlamLinearSolver>();
    // auto linearSolver = g2o::make_unique<CholmodLinearSolver>();
    auto linearSolver = g2o::make_unique<CSparseLinearSolver>();
    linearSolver->setBlockOrdering(false);
    g2o::OptimizationAlgorithmGaussNewton* solver = new g2o::OptimizationAlgorithmGaussNewton(
        g2o::make_unique<SlamBlockSolver>(std::move(linearSolver)));
    // g2o::OptimizationAlgorithmLevenberg* solver = new g2o::OptimizationAlgorithmLevenberg(
    //     g2o::make_unique<SlamBlockSolver>(std::move(linearSolver)));

    optimizer.setAlgorithm(solver);

    // ********************************************************
    // Building graph: 
    //      - Adding poses/vertices
    //      - Adding odometry measurements
    if(verb)
        std::cout << "Building factor graph..." << std::endl;

    // Index that keeps track of the gps measurements
    size_t idx_gps = 0;
    for(int k = 0; k < K; k++){
        g2o::SE2::VertexSE2* robot = new g2o::SE2::VertexSE2;
        robot->setId(k);
        // Set estimate from the L-InEKF
        robot->setEstimate( X_kf[k]);
        // robot->setEstimate( Pose::Identity());
        optimizer.addVertex(robot);

        // Add prior measurement
        if( k == 0){
            g2o::SE2::UEdgeSE2Prior* prior = new g2o::SE2::UEdgeSE2Prior;
            prior->vertices()[0] = optimizer.vertex( k);                        
            EPose X_prior_eig = meas_prior.mean();
            Pose  X_prior = Pose( X_prior_eig(0, 2), X_prior_eig( 1, 2), std::atan2( X_prior_eig(1, 0), X_prior_eig( 0, 0)));
            prior->setMeasurement( X_prior);
            prior->setInformation( meas_prior.cov().inverse());
        }
        if(k > 0){
            int km1 = k - 1;
            // Add odometry edges
            g2o::SE2::BEdgeSE2SE2* odom = new g2o::SE2::BEdgeSE2SE2;
            // odom->vertices()[0] = optimizer.vertex( k - 1);
            // odom->vertices()[1] = optimizer.vertex( k);
            odom->setVertex( 0, optimizer.vertex( k - 1));
            odom->setVertex( 1, optimizer.vertex( k));
            // Compute the measurement vector
            double dt_km1 = dt_func(km1);
            odom->setDt( dt_km1);
            odom->setMeasurement( 
                g2o::Vector3( 
                    meas_vel [km1].mean()(0),
                    meas_vel [km1].mean()(1),
                    meas_gyro[km1].mean()(0)                    
                )
            );
            // Compute process noise covariance
            CovQ Q_km1  = CovQ::Zero();
            Q_km1.block< dof_vel, dof_vel>(0, 0)   = meas_vel [km1].cov();
            Q_km1.block< dof_gyro, dof_gyro>(2, 2) = meas_gyro[km1].cov();
            // Jacobian of process model w.r.t. process noise w_km1
            JacF_wkm1 jac_F_wkm1 = dt_km1 * JacF_wkm1::Identity();
            odom->setInformation( (jac_F_wkm1 * Q_km1 * jac_F_wkm1.transpose()).inverse());
            
            // Add edge to graph
            optimizer.addEdge( odom);
        }

        // Add GPS edges
        if( idx_gps < meas_gps.size() && (meas_gps[idx_gps].time() <= X_kf_rv[k].time())){
            // Create GPS unary edge
            g2o::SE2::UEdgeSE2Gps* y_gps_k = new g2o::SE2::UEdgeSE2Gps;
            // y_gps_k->vertices()[0] = optimizer.vertex( k);
            y_gps_k->setVertex( 0, optimizer.vertex( k));
            // Set the measurement
            y_gps_k->setMeasurement(
                    g2o::Vector2(
                        meas_gps[ idx_gps].mean()
                    )
                );
            // TODO: Compute the Jacobian of the gps-error w.r.t. state
            // Compute Jacobian of the gps-error w.r.t. measurement noise.
            JacYgps_nk jac_Y_gps_n = X_kf[k].rotation().transpose();
            // Compute GPS noise covariance
            CovGps R_k = meas_gps[idx_gps].cov();
            // Set information matrix
            y_gps_k->setInformation( (jac_Y_gps_n * R_k * jac_Y_gps_n.transpose()).inverse());

            // Add edge to graph
            optimizer.addEdge( y_gps_k);

            // Increment GPS measurement counter
            idx_gps++;
        }
    }

    // dump initial state to the disk
    optimizer.save("tutorial_before.g2o");

    // ********************************************************
    // Solving the optimization problem
    if(verb)
        std::cout << "Setting up optimization problem" << std::endl;

#ifndef NDEBUG
    if(verb)
        std::cout << "\tFixing first pose" << std::endl;
#endif
    // fix the first robot pose to account for gauge freedom
    g2o::SE2::VertexSE2* firstRobotPose = dynamic_cast<g2o::SE2::VertexSE2*>(optimizer.vertex(0));
    firstRobotPose->setFixed( config["g2o"]["fix_x0"].as<bool>());
    // Note: if fixing pose, then use a ground truth prior!
    if(verb)
        std::cout << "First pose fixed: " << optimizer.vertex(0)->fixed() << std::endl;

    // Set verbosity
    optimizer.setVerbose( verb);

    if(verb)
        std::cout << "Optimizing" << std::endl;
    optimizer.initializeOptimization();
    // optimizer.computeInitialGuess();
    optimizer.optimize( config["g2o"]["max_iterations"].as<int>());
    // optimizer.computeBatchStatistics();
    if(verb)
        std::cout << "Done" << std::endl;


    // ********************************************************
    // Computing marginals    
    if(verb)
        std::cout << "Computing marginals" << std::endl;
    // Compute marginals
    auto vertices = optimizer.activeVertices();
    
    // ********************************************************
    // Export to random variable vector    
    std::vector< PoseEstimate> X_batch_rv( K);
    // Progress bar
    // Show an example of an "empty" progress bar
    if(verb){
        std::cout << "|" << std::string( 61, ' ') << "|" << std::endl;
        std::cout << "|" << std::flush;
    }
    for( int k = 0; k < K; k++){
        // Output iteration
        if( verb && k % (K/60) == 0){
            // std::cout << "\tMarginals: " << k << " of " << K << std::endl;
            std::cout << "=" << std::flush;
        }
        // Set time
        X_batch_rv[k].setTime( X_kf_rv[k].time());

        // Set mean
        g2o::SE2::VertexSE2* p_X_k = dynamic_cast<g2o::SE2::VertexSE2*>(optimizer.vertex(k));
        X_batch_rv[k].setMean( p_X_k->estimate().transform());       
        
        if( firstRobotPose->fixed()){
            if( k == 0){
                X_batch_rv[0].setCov( 
                    0 * X_kf_rv[0].cov()
                ); 
            }else{
                // Set covariance        
                g2o::SparseBlockMatrixX sp_cov_X_k;
                // Vertex container containing 1 vertex
                g2o::OptimizableGraph::VertexContainer v_k( 1); 
                v_k[0] = vertices[k];
                optimizer.computeMarginals( sp_cov_X_k, v_k);
                if(sp_cov_X_k.rows() > 0){
                    X_batch_rv[k].setCov( 
                        CovPosThetaToCovThetaPos( sp_cov_X_k.block( k - 1, k - 1)->eval())
                    );
                }
            }
        }else{
            // Set covariance        
            g2o::SparseBlockMatrixX sp_cov_X_k;
            // Vertex container containing 1 vertex
            g2o::OptimizableGraph::VertexContainer v_k( 1); 
            // Set covariance        
            v_k[0] = vertices[k];
            optimizer.computeMarginals( sp_cov_X_k, v_k);
            if(sp_cov_X_k.rows() > 0){
                optimizer.computeMarginals( sp_cov_X_k, v_k);
                
                X_batch_rv[k].setCov( 
                    CovPosThetaToCovThetaPos( sp_cov_X_k.block( k, k)->eval())
                );
            }
        }
    }
    // End of progress bar
    if(verb)
        std::cout << "|" << std::endl;

    // Exporting KF estimates
    if(verb)
        std::cout << "Exporting L-InEKF estimate to " << config["filename_kf"].as<std::string>() << std::endl;
    RV::IO::write( X_kf_rv, config["filename_kf"].as<std::string>(), "X");

    // Exporting batch estimates
    if(verb)
        std::cout << "Exporting batch estimate to " << config["filename_batch"].as<std::string>() << std::endl;
    RV::IO::write( X_batch_rv, config["filename_batch"].as<std::string>(), "X");

    // Free graph memory
    optimizer.clear();
}
