function [X_hat, P_hat] = initLinekf(struct_prior, struct_vel, struct_gyro, struct_gps, t_sim)
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
    % ----------------------------------
    % Output:
    %   X_hat   :   [ 3 x 3 x K] double
    %       3D array of the L-InEKF estiamtes
    %   P_hat   :   [ 3 x 3 x K]
    %       Covariance on these estimates
    %
    % ----------------------------------
    %   Note that this requires the internal MLG package.
    %
    %   Amro Al-Baali
    %   08-May-2021
    % ----------------------------------
    
    
    % Get the number of poses
    K = length( t_sim);
    
    % Measurements
    %   Creat the odometry array: u_k = [ gyro_k; vel_k];
    u_arr = [ struct_gyro.mean; struct_vel.mean];
    %   GPS measurements
    y_gps = struct_gps.mean;
    t_gps = struct_gps.time;
    n_gps = length( t_gps);
    
    
    % Create arrays
    X_hat = zeros( 3, 3, K);
    P_hat = zeros( 3, 3, K);
    
    % Initialize with the prior
    X_hat( :, :, 1) = struct_prior.mean;
    P_hat( :, :, 1) = struct_prior.cov;
    
    % Latest GPS measurement
    idx_gps_j = 1;
    for kk = 2 : K
        % Get parameters at kk - 1
        dt_km1 = t_sim( kk) - t_sim( kk - 1);
        u_km1 = u_arr( :, kk - 1);
        X_km1 = X_hat( :, :, kk - 1);
        P_hat_km1 = P_hat( :, :, kk - 1);
        
        % Process noise (the dt_km1^2 is from L * () * L', where L = dt * eye.
        Q_km1 = dt_km1^2 * blkdiag( struct_gyro.cov( :, :, kk - 1), ...
                                    struct_vel.cov( :, :, kk - 1));
        % Jacobian of process model w.r.t. previous state
        A_km1 = SE2.adjoint( se2alg.expMap( - dt_km1 * u_km1));
        % Predict
        X_check_k = X_km1 * se2alg.expMap( dt_km1 * u_km1);
        P_check_k = A_km1 * P_hat_km1 * A_km1' + Q_km1;
        
        % Check if there's a correction
        if idx_gps_j <= n_gps && ( t_gps( idx_gps_j) <= t_sim( kk))
            % Measurement noise covariance
            R_k = struct_gps.cov( :, :, idx_gps_j);
            % Measurement
            y_k = y_gps( :, idx_gps_j);
            % Left-invariant innovation
            z_k = X_check_k( 1 : 2, 1 : 2)' * ( y_k - X_check_k( 1 : 2, 3));
            % Jacobians of measurement model 
            %   w.r.t. state
            H_k = - [ zeros( 2, 1), eye( 2)];
            %   w.r.t. measurement noise = C_ab_check_k'
            M_k = X_check_k( 1 : 2, 1 : 2)';
            
            % Compute Kalman gain
            K_k = P_check_k * H_k' / ( H_k * P_check_k * H_k' + M_k * R_k * M_k');
            
            % Correct
            X_hat_k = X_check_k * se2alg.expMap( - ( K_k * z_k));
            P_hat_k = (eye(3) - K_k * H_k) * P_check_k * ...
                    (eye(3) - K_k * H_k)'  + K_k * M_k * R_k * M_k' * K_k';
            
            % Ensure symmetry of the covariance
            P_hat_k = (1/2) * ( P_hat_k + P_hat_k');
            
            % Increment gps index
            idx_gps_j = idx_gps_j + 1;
        else
            X_hat_k = X_check_k;
            P_hat_k = P_check_k;
        end
        
        % Add corrections
        X_hat( :, :, kk) = X_hat_k;
        P_hat( :, :, kk) = P_hat_k;
    end    
end