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
% clear all;
% close all;

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
% Load noisy measurements (.mat file)
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

%% Run L-InEKF
[ X_kf, P_kf] = Initialization.initLinekf( data_struct.meas.prior, ...
                                    data_struct.meas.velocity, ...
                                    data_struct.meas.gyro, ...
                                    data_struct.meas.gps, ...
                                    t_sim);
%% Batch initialization
fprintf("Initializing states using '%s'\n", config_yml.init_method);
tic();
switch lower(config_yml.init_method)
    case 'odom'
        X_initial = Initialization.initOdom( X_prior, meas_vel, meas_gyro, t_sim);
    case 'l-inekf'   
        X_initial = X_kf;
        
end
disp('Done');
toc();

%% Set up the batch problem
[ X_batch, infm_batch ] = batchOptimization.batchOptimizationSE2(data_struct.meas.prior, ...
                                    data_struct.meas.velocity, ...
                                    data_struct.meas.gyro, ...
                                    data_struct.meas.gps, ...
                                    t_sim, ...
                                    X_initial, ...
                                    config_yml.optim_params);


%% Analysis
K = size( X_batch, 3);

% Get ground truth
X_gt = load( config_yml.filename_gt).X_poses;

% Compute errors
xi_kf_arr = cell2mat( arrayfun( @(kk) SE2.Log( X_gt(:, :, kk) \ X_kf(:, :, kk)), ...
    1 : K, 'UniformOutput', false));
xi_batch_arr = cell2mat( arrayfun( @(kk) SE2.Log( X_gt(:, :, kk) \ X_batch(:, :, kk)), ...
    1 : K, 'UniformOutput', false));

% Get batch covariances
P_batch = getBlockDiagonals( inv( infm_batch), 3);

%% Plots
% Time
t_sim = data_struct.sim.time;
% Colors
col_kf_err = matlabColors( 'orange');
col_kf_var = col_kf_err; %matlabColors( 'orange');

col_batch_err = matlabColors( 'blue');
col_batch_var = col_batch_err; %matlabColors( 'orange');

% Error plots
figure;
for lv1 = 1 : 3
    subplot( 3, 1, lv1);
    hold all; grid on;
    % Filter
    plot( t_sim, xi_kf_arr( lv1, :), 'LineWidth', 1.5, 'Color', col_kf_err, ...
        'DisplayName', 'Filter');
    plot( t_sim, 3 * sqrt( squeeze( P_kf( lv1, lv1, :))), '-.', 'LineWidth', 1.5, ...
        'Color', col_kf_var, 'HandleVisibility', 'off');
    plot( t_sim, -3 * sqrt( squeeze( P_kf( lv1, lv1, :))), '-.', 'LineWidth', 1.5, ...
        'Color', col_kf_var, 'HandleVisibility', 'off');
    
    % Batch
    plot( t_sim, xi_batch_arr( lv1, :), 'LineWidth', 1.5, 'Color', col_batch_err, ...
        'DisplayName', 'Batch');
    plot( t_sim, 3 * sqrt( squeeze( P_batch( lv1, lv1, :))), '-.', 'LineWidth', 1.5, ...
        'Color', col_batch_var, 'HandleVisibility', 'off');
    plot( t_sim, -3 * sqrt( squeeze( P_batch( lv1, lv1, :))), '-.', 'LineWidth', 1.5, ...
        'Color', col_batch_var, 'HandleVisibility', 'off');
    
    ylabel(sprintf('$\\delta\\xi_{%i}$', lv1), 'Interpreter', 'latex', 'FontSize', 14);
    if lv1 == 1
        legend('Interpreter', 'latex', 'FontSize', 14);
    end
end
xlabel('$t_{k}$ [s]', 'Interpreter', 'latex', 'FontSize', 14);

%% Plot trajectory
% Get states
X_gt_states( K)     = StateSE2();
X_kf_states( K)     = StateSE2(); 
X_batch_states( K) = StateSE2();

for kk = 1 : K
    % Ground truth
    X_gt_states( kk).state = X_gt( :, :, kk);
    X_gt_states( kk).time  = t_sim( kk);
    
    % Filter estimates
    X_kf_states( kk).state = X_kf( :, :, kk);
    X_kf_states( kk).time  = t_sim( kk);
    
    % Batch estimates
    X_batch_states( kk).state = X_batch( :, :, kk);
    X_batch_states( kk).time  = t_sim (kk);
end

col_kf    = matlabColors( 'orange');
col_batch = matlabColors( 'blue');

figure; 
plotMlgPose( X_gt_states, '-', matlabColors('grey'));
hold on;
plotMlgPose( X_kf_states, '-.', col_kf);
plotMlgPose( X_batch_states, '-.', col_batch);
legend({'Ground truth', 'L-InEKF', 'Batch'}, 'Interpreter', 'latex', 'FontSize', 14);
