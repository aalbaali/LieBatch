%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   L-InEKF SE(2) batch SLAM using the following measurements:
%       - linear velocity,
%       - angular velocity (i.e., rate gyro), and
%       - GPS.
%   The noisy data should be provided, and the output is the state estimates and
%   covariances. 
%
%   No plots are provded, thus not ground-truth is needed.
%
%   Amro Al-Baali
%   08-May-2021
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;
close all;

%% Settings
% Add project paths
disp('Make sure the program is running from the SE2 directory (not the root directory)');
addprojectpaths();


% Load paths from YAML config file (YAML reader is loaded in the extern files)
% The configuration file should contain the batch initialization method. 
% The options so far are
%   - odom (odometry)
%   - L-InEKF (Left-invariant EKF)
config_yml = YAML.read( 'config.yml');

%% Load data
data_struct = load( config_yml.filename_data).data_struct;

% Simulation
t_sim = data_struct.sim.time;

% Measurements
%   Prior
X_prior = data_struct.meas.prior.mean;
P_prior = data_struct.meas.prior.cov;
%   Linear velocity
meas_vel = data_struct.meas.velocity.mean;
cov_vel  = data_struct.meas.velocity.cov;
%   Angular velocity
meas_gyro = data_struct.meas.gyro.mean;
cov_gyro  = data_struct.meas.gyro.cov;
%   GPS
meas_gps = data_struct.meas.gps.mean;
cov_gps  = data_struct.meas.gps.cov;
%% Initialize
fprintf("Initializing states using '%s'\n", config_yml.init_method);
tic();
switch lower(config_yml.init_method)
    case 'odom'
        X_initial = initOdom( X_prior, meas_vel, meas_gyro, t_sim);
    case 'l-inekf'
end
disp('Done');
toc();