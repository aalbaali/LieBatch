function batchOptimizationSE2(struct_prior, struct_vel, struct_gyro, struct_gps, ...
    t_sim, X_initial, params)
    %INITLINEKF(struct_prior, struct_vel, struct_gyro, struct_gps, t_sim)
    %generates left-invariant EKF solution.
    %
    % Input:
    %   struct_prior    :   struct
    %       SE(2) prior at time t_sim(1). 
    %       Contains    
    %           meas        :   [ 3 x 3] double
    %           cov         :   [ 3 x 3] SPD matrix    
    %   struct_vel      :   struct
    %       Linear velocities, where K is the number of poses.
    %       Contains
    %           meas        :   [ 2 x K - 1] double
    %           cov         :   [ 2 x 2 x K - 1] SPD matrices
    %           time        :   [ 1 x K - 1] measurements time
    %   struct_gyro      :   struct
    %       Angular velocities, where K is the number of poses.
    %       Contains
    %           meas        :   [ 1 x K - 1] double
    %           cov         :   [ 1 x 1 x K - 1] SPD matrices
    %           time        :   [ 1 x K - 1] measurements time
    %   struct_gps      :   struct
    %       GPS measurements, where n_gps is the number of GPS meas.
    %       Contains
    %           meas        :   [ 2 x n_gps] double
    %           cov         :   [ 2 x 2 x n_gps] SPD matrices    
    %           time        :   [ 1 x n_gps] measurements time
    %   t_sim       :   [ 1 x K] double
    %       Simulation time steps
    %   X_initial   :   [ 3 x 3 x K]
    %       Optimization initial guess.
    %   params  :   struct
    %       Optimization parameters (loaded from the config.yml file).
    %       
    % ----------------------------------
    % Output:
    %   X_hat   :   [ 3 x 3 x K] double
    %       3D array of the L-InEKF estiamtes    
    %   P_hat   :   [ 3 x 3 x K]
    %       Marginalized covariance on these estimates
    %   P_hat_joint : [ 3K x 3K]
    %       Joint covariance on all X_hat
    %
    % ----------------------------------
    %   Note that this requires the internal MLG package.
    %
    %   Amro Al-Baali
    %   08-May-2021
    % ----------------------------------
    
    % Flag to indicate whether the optimization succeded or not
    successful = false;
    
    % Number of poses
    K = length( t_sim);
    
    
    % Measurements
    %   Prior
    X_prior = struct_prior.mean;
    %   Odometry
    %   Creat the odometry array: u_k = [ gyro_k; vel_k];
    u_arr = [ struct_gyro.mean; struct_vel.mean];
    %   GPS measurements
    y_gps   = struct_gps.mean;
    t_gps   = struct_gps.time;
    n_gps   = length( t_gps);
    idx_gps = ceil( t_gps / ( t_sim( 2) - t_sim( 1))); % Assuming frequency is constant
    
    dd =errorFunction( X_initial, struct_prior, struct_vel, struct_gyro, struct_gps, t_sim);
end


function [ err_val, err_jac] = errorFunction( X, struct_prior, struct_vel, struct_gyro, struct_gps, t_sim)
    % Batch error function
    
    % Number of poses
    K = length( t_sim);
    % Number of states per pose (constant value)
    n_x = 3;
    
    % Measurements
    %   Prior
    X_prior = struct_prior.mean;
    %   Odometry
    %   Creat the odometry array: u_k = [ gyro_k; vel_k];
    u_arr = [ struct_gyro.mean; struct_vel.mean];
    %   GPS measurements
    y_gps   = struct_gps.mean;
    t_gps   = struct_gps.time;
    n_gps   = size( y_gps, 1);
    K_gps   = length( t_gps);
    idx_gps = ceil( t_gps / ( t_sim( 2) - t_sim( 1))); % Assuming frequency is constant
    
    % Set up the arrays and matrices
    %   Prior
    err_prior   = SE2.Log( X( :, :, 1) \ X_prior);
    jac_prior_x = kron( sparse( 1, 1, 1, 1, K), speye( n_x));
    %   Odometry 
    %       Error array
    err_odom = nan( n_x, K - 1);
    %       Jacobians w.r.t. states
    jac_odom_x = sparse( [], [], [], 0, n_x * K);


    %   GPS
    %       Error array
    err_gps = nan( n_gps, K_gps);
    %       Jacobian of GPS err w.r.t. states
    jac_gps_x = sparse( [], [], [], 0, n_x * K);
    
    % Latest GPS measurement
    idx_gps_j = 1;
    % Let's build the matrices!
    for kk = 1 : K
        % Odometry
        if kk >= 2
            % Sampling period
            dt_km1 = t_sim( kk) - t_sim( kk - 1);
            % Left-invariant error
            %   X_check_k = F_km1( X_km1, u_km1);
            X_check_k = X( :, :, kk - 1) * se2alg.expMap( dt_km1 * u_arr( :, kk - 1));
            err_odom( :, kk - 1) = SE2.Log( X( :, :, kk) \ X_check_k);
            % Jacobian of the current error function w.r.t. the states
            jac_odom_x_k = kron( sparse( 1, kk, 1, 1, K), speye( n_x)) ...
                + kron( sparse( 1, kk - 1, 1, 1, K), ...
                    - SE2.adjoint( se2alg.expMap( -dt_km1 * u_arr( :, kk - 1))));
            %   Jacobians
            jac_odom_x = [ jac_odom_x;
                           jac_odom_x_k];
                        
        end
        
        % Check for exteroceptive meas.
        if idx_gps_j <= K_gps && ( t_gps( idx_gps_j) <= t_sim( kk))
            % Measurement
            y_k = y_gps( :, idx_gps_j);
            
            % Error
            err_gps( :, idx_gps_j) = X( 1 : 2, 1 : 2, kk)' * ( y_k - X( 1 : 2, 3, kk));
            % Jacobians of measurement model 
            %   w.r.t. state
            H_k = - [ zeros( 2, 1), eye( 2)];
            jac_gps_x = [ jac_gps_x;
                    kron( sparse( 1, kk, 1, 1, K), H_k)];
            
            % Increment gps index
            idx_gps_j = idx_gps_j + 1;
        end
    end
    
    % Augment the error arrays and Jacobians
    err_val = [ err_prior; reshape( err_odom, [], 1); reshape( err_gps, [], 1)];
    err_jac = [ jac_prior_x; jac_odom_x; jac_gps_x];    
end