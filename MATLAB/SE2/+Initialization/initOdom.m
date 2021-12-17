function [ X_check, P_check] = initOdom( meas_prior, meas_vel, meas_gyro, t_sim)
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
  %   X_check   :   [ 3 x 3 x K] double
  %       3D array of the dead-reckoning solution
  %   P_check   :   [ 3 x 3 x K] SDP covariance matrix
  %       Dead-reckoning covariance matrix
  %
  % ----------------------------------
  %   Note that this requires the internal MLG package.
  %
  %   Amro Al-Baali
  %   08-May-2021
  % ----------------------------------

  % Get the number of poses
  K = length( t_sim);

  % Creat the odometry array: u_k = [ gyro_k; vel_k];
  u_arr = [ meas_gyro.mean; meas_vel.mean];

  % Create array
  X_check = zeros( 3, 3, K);
  P_check = zeros( 3, 3, K);

  % Initialize with the prior
  X_check( :, :, 1) = meas_prior.mean;
  P_check( :, :, 1) = meas_prior.cov;

  % Dead-reckoning
  for kk = 2 : K
    % Sampling period
    dt_km1 = t_sim( kk) - t_sim( kk - 1);
    % Xi
    Xi_km1 = se2alg.expMap( dt_km1 * u_arr( :, kk - 1));
    % Propagate the mean
    X_check( :, :, kk) = X_check( :, :, kk - 1) * Xi_km1;

    % Error
    e_k = SE2.Log( X_check( :, :, kk) \ X_check( :, :, kk - 1) * Xi_km1);
    % Jacobians
    J_xkm1 = -SE2.computeJRightInv( e_k) * SE2.adjoint( SE2.inverse( Xi_km1));
    L_km1  = -SE2.computeJRightInv( e_k) * dt_km1;
    % Process noise covariance
    Q_km1  = blkdiag( meas_gyro.cov( :, :, kk - 1), meas_vel.cov( :, :, kk - 1));
    % Propagate covariances
    P_check( :, :, kk) = J_xkm1 * P_check( :, :, kk - 1) * J_xkm1' + L_km1 * Q_km1 * L_km1';
  end
end