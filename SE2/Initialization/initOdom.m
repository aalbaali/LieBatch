function X_hat = initOdom( meas_prior, meas_vel, meas_gyro, t_sim)
    %INITODOM(meas_prior, meas_vel, meas_gyro, t_sim) generates dead-reckoning
    %solution using odometry only.
    %
    % Input:
    %   meas_prior  :   [ 3 x 3] double 
    %       SE(2) prior at time t_sim(1)
    %   meas_vel    :   [ 2 x K - 1] double
    %       Linear velocities, where K is the number of poses
    %   meas_gyro   :   [ 1 x K - 1] double
    %       Angular velocities, where K is the number of poses
    %   t_sim       :   [ 1 x K] double
    %       Simulation time steps
    % ----------------------------------
    % Output:
    %   X_hat   :   [ 3 x 3 x K] double
    %       3D array of the dead-reckoning solution
    %
    % ----------------------------------
    %   Note that this requires the internal MLG package.
    %
    %   Amro Al-Baali
    %   08-May-2021
    
    % Get the number of poses
    K = length( t_sim);
    
    % Creat the odometry array: u_k = [ gyro_k; vel_k];
    u_arr = [ meas_gyro; meas_vel];
    
    % Create array
    X_hat = zeros( 3, 3, K);
    
    % Initialize with the prior
    X_hat( :, :, 1) = meas_prior;
    
    % Dead-reckoning
    for kk = 2 : K
        % Sampling period
        dt_km1 = t_sim( kk) - t_sim( kk - 1);
        X_hat( :, :, kk) = X_hat( :, :, kk - 1) ...
            * se2alg.expMap( dt_km1 * u_arr( :, kk - 1));
    end
end