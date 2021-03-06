%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   L-InEKF SE(2) batch SLAM using the following measurements:
%       - linear velocity,
%       - angular velocity (i.e., rate gyro), and
%       - GPS.
%   The noisy data should be provided, and the output is the state estimates and
%   covariances. 
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
% config_yml = YAML.read( 'config.yml');
config_yml = ReadYaml('config.yml');

data_dir_out = 'G:\My Drive\Professional\Code_base\Local\MATLAB\Reseach_codebase\McGill_2G\delayed_batch\Data_lessnoisy_2\Processed_full_data';
% data_dir_out = 'G:\My Drive\Professional\Code_base\Local\MATLAB\Research_codebase\McGill_2G\delayed_batch\Data_lessnoisy\Processed_full_data_BatchPrior_onFirstLC';

%% Load data
% Load noisy measurements (.mat file)
data_struct = load( config_yml.filename_data).data_struct;

% Simulation
t_sim = data_struct.sim.time;

% Batch parameters
batch_params = config_yml.batch_params;

% Measurements
%   Prior
prior_struct = data_struct.meas.prior;
%   Linear velocity
meas_vel_struct = data_struct.meas.velocity;
%   Angular velocity
meas_gyro_struct = data_struct.meas.gyro;
%   GPS
meas_gps_struct = data_struct.meas.gps;

%   LC (if it exists)
if batch_params.include_lc && isfield( data_struct.meas, 'lc')
  include_lc = batch_params.include_lc;
  % Measurements
  lcs_struct = data_struct.meas.lc;
else
  % Don't include LC if the measurements don't exist
  warning('No LC is included');
  include_lc = false;
  lcs_struct = struct( 'mean', [], 'cov', [], 'idx', [], 'time', []);
end
%% Run L-InEKF
[ X_kf, P_kf] = Initialization.initLinekf( data_struct.meas.prior, ...
                                    data_struct.meas.velocity, ...
                                    data_struct.meas.gyro, ...
                                    meas_gps_struct, ...
                                    t_sim);
%% Batch initialization
fprintf("Initializing states using '%s'\n", config_yml.init_method);
tic();
switch lower(config_yml.init_method)
    case 'odom'
        [ X_initial, P_odom] = Initialization.initOdom( prior_struct, meas_vel_struct, meas_gyro_struct, t_sim);
    case 'l-inekf'   
        X_initial = X_kf;
end
disp('Done');
toc();

%% Set up the batch problem
[ X_batch, infm_batch ] = batchOptimization.batchOptimizationSE2(data_struct.meas.prior, ...
                                    data_struct.meas.velocity, ...
                                    data_struct.meas.gyro, ...
                                    meas_gps_struct, ...
                                    lcs_struct, ...
                                    t_sim, ...
                                    X_initial, ...
                                    config_yml.batch_params, ...
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
X_initial_states    = StateSE2();
X_gt_states( K)     = StateSE2();
X_kf_states( K)     = StateSE2(); 
X_batch_states( K)  = StateSE2();

for kk = 1 : K
% Initial estimate
X_initial_states( kk).state = X_initial(:, :, kk);
X_initial_states( kk).time  = t_sim( kk);

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

% Save batch
if strcmpi( config_yml.init_method, 'odom')
  X_odom_states = X_initial_states;
else
  warning('X_initial is the IEKF initialization');
end
save( fullfile( data_dir_out, 'full_odom'), 'X_odom_states', 'P_odom')
save( fullfile( data_dir_out, 'full_kf'), 'X_kf_states', 'P_kf')
save( fullfile( data_dir_out, 'full_batch'), 'X_batch_states', 'P_batch')
save( fullfile( data_dir_out, 'full_gt'), 'X_gt_states');
fprintf("Saved files to\n\t%s\n", data_dir_out);