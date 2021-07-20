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

% data_dir_out = 'G:\My Drive\Professional\Code_base\Local\MATLAB\Reseach_codebase\McGill_2G\delayed_batch\Data_lessnoisy_2\Processed_full_data';
data_dir_out = 'G:\My Drive\Professional\McGill\Masters\Data\2G\Simulation\figure_8_trajectory_lessnoisy_mct\processed_full';

%%%%%%%%%%%%%%%%
warning('TEMPORARY: Using data directory and ignoring data struct from the config file');
data_dir = 'G:\My Drive\Professional\McGill\Masters\Data\2G\Simulation\figure_8_trajectory_lessnoisy_mct\noisy_data';
% Get the noisy data files
noisy_files = dir( fullfile( data_dir, 'noisy_data*.mat'));
%%%%%%%%%%%%%%%%
%% Load data
% % Load noisy measurements (.mat file)
% data_struct = load( config_yml.filename_data).data_struct;

num_trials = length( noisy_files);
warning('TEMPORARY: Using data directory and ignoring data struct from the config file');
% MCT timer
tic_mct  = tic();
for lv_trial = 1 : 2 % num_trials
  fprintf('Starting trial\t%i\tof\t%i\n', lv_trial, num_trials);
  
  fprintf('\tLoading data...\n');
  % Filename
  filename_noisy = noisy_files( lv_trial).name;
  % Load data struct
  data_struct = load( fullfile( data_dir, filename_noisy)).data_struct;

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
  fprintf('\tRunning InEKF data...\n');
  [ X_kf, P_kf] = Initialization.initLinekf( data_struct.meas.prior, ...
                                      data_struct.meas.velocity, ...
                                      data_struct.meas.gyro, ...
                                      meas_gps_struct, ...
                                      t_sim);
                                    
  %% Batch initialization
  fprintf("\tInitializing states using '%s'\n", config_yml.init_method);
%   tic();
  switch lower(config_yml.init_method)
      case 'odom'
          [ X_initial, P_odom] = Initialization.initOdom( prior_struct, meas_vel_struct, meas_gyro_struct, t_sim);
      case 'l-inekf'   
          X_initial = X_kf;
  end
  disp('\tDone');
%   toc();

  %% Set up the batch problem
  fprintf('\tRunning batch...\n');
  fprintf('\t');
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

  % Get batch covariances
  P_batch = getBlockDiagonals( inv( infm_batch), 3);

  %% Plots
  % Time
  t_sim = data_struct.sim.time;

  %% Trajectories
  fprintf('\tGenerating trajectories...\n');
  % Get states
  X_initial_states    = StateSE2();
  if lv_trial == 1
    X_gt_states( K)     = StateSE2();
  end
  X_kf_states( K)     = StateSE2(); 
  X_batch_states( K)  = StateSE2();

  for kk = 1 : K
    % Initial estimate
    X_initial_states( kk).state = X_initial(:, :, kk);
    X_initial_states( kk).time  = t_sim( kk);

    if lv_trial == 1
      % Ground truth
      X_gt_states( kk).state = X_gt( :, :, kk);
      X_gt_states( kk).time  = t_sim( kk);
    end

    % Filter estimates
    X_kf_states( kk).state = X_kf( :, :, kk);
    X_kf_states( kk).time  = t_sim( kk);

    % Batch estimates
    X_batch_states( kk).state = X_batch( :, :, kk);
    X_batch_states( kk).time  = t_sim (kk);
  end

  
  %% Save batch
  fprintf('\tSaving data\n');
  if strcmpi( config_yml.init_method, 'odom')
    X_odom_states = X_initial_states;
  else
    warning('X_initial is the IEKF initialization');
  end
  
  save( fullfile( data_dir_out, sprintf('full_odom_%02i', lv_trial)), ...
    'X_odom_states', 'P_odom')
  save( fullfile( data_dir_out, sprintf('full_kf_%02i', lv_trial)), 'X_kf_states', 'P_kf')
  save( fullfile( data_dir_out, sprintf('full_batch_%02i', lv_trial)), 'X_batch_states', 'P_batch')  
end
toc( tic_mct);
save( fullfile( data_dir_out, 'full_gt'), 'X_gt_states');
fprintf("Saved files to\n\t%s\n", data_dir_out);