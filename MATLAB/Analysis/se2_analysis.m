%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Importing KF estimates and plotting error plots.
%
%   Amro Al Baali
%   15-Mar-2021
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% clear all;
% close all;

%% Load data
% Ground truth file name
file_name_gt = "\\wsl$\Ubuntu-20.04\home\aa\Documents\Data\Data_generator\SE2\gt_states.txt";
% IEKF state estimates file name
file_name_kf    = "\\wsl$\Ubuntu-20.04\home\aa\Documents\Data\Data_generator\SE2\est_xhat.txt";
% Batch state estimates file name
file_name_batch = "\\wsl$\Ubuntu-20.04\home\aa\Documents\Data\Data_generator\SE2\est_xhat_batch.txt";

% Load data
%   Specify the number of columns since the ground truth doesn't include
%   covariances
struct_gt    = importTextFile( file_name_gt);
struct_kf    = importTextFile( file_name_kf);
struct_batch = importTextFile( file_name_batch);

% Plotting range
idx_range = 1 : min( length(struct_gt.time), length( struct_batch.time));

plt_kf    = true;
plt_batch = true;

%% Normalize SE(2) objects
% Re-normalize poses lambda function
renormalize = @(X_k) se2alg.expMap( se2alg.vee( SE2.logMap( X_k)));
X_kf    = struct_kf.values3d;
X_batch = struct_batch.values3d;
X_gt    = struct_gt.values3d;
% Create array of StateSE2 objects
K = size( X_batch, 3);
X_kf_states( K)    = StateSE2();
X_batch_states( K) = StateSE2();
X_gt_states( K)    = StateSE2();

% Normalize each pose
for kk = 1 : length( idx_range)
    % Ground-truth
    C_k = X_gt( 1 : 2, 1 : 2, kk);
    C_k = SO2.synthesize( atan2( C_k(2, 1), C_k(1, 1)));
    X_gt( 1 : 2, 1 : 2, kk) = C_k;
    X_gt_states( kk).state = X_gt( :, :, kk);
    X_gt_states( kk).time  = struct_gt.time( kk);
    
    % Filter
    C_k = X_kf( 1 : 2, 1 : 2, kk);
    C_k = SO2.synthesize( atan2( C_k(2, 1), C_k(1, 1)));
    X_kf( 1 : 2, 1 : 2, kk) = C_k;
    % StateSE2 objects
    X_kf_states( kk).state = X_kf( :, :, kk);
    X_kf_states( kk).time  = struct_kf.time( kk);
    
    % Batch
    C_k = X_batch( 1 : 2, 1 : 2, kk);
    C_k = SO2.synthesize( atan2( C_k(2, 1), C_k(1, 1)));
    X_batch( 1 : 2, 1 : 2, kk) = C_k;
    % StateSE2 objects
    X_batch_states( kk).state = X_batch( :, :, kk);
    X_batch_states( kk).time  = struct_batch.time( kk);
end

%% Errors
% Lambda function to compute error
SE2Ominus = @(X_k1, X_k2) se2alg.vee(SE2.logMap( X_k1 \ X_k2));
% Compute errors
%   Filter
dxi_kf    = cell2mat( arrayfun(@(kk) SE2Ominus( X_gt(:,:,kk), ...
    X_kf(:,:,kk)), idx_range, 'UniformOutput', false));
%   Batch
dxi_batch = cell2mat( arrayfun(@(kk) SE2Ominus( X_gt(:,:,kk), ...
    X_batch(:,:,kk)), idx_range, 'UniformOutput', false));
% Get covariance matrices
time = struct_gt.time;

%% Plots
figure; 

% Covariances
%   Filter
P_kf    = struct_kf.cov;
%   Batch 
P_batch = struct_batch.cov;

% % Convert from local to global covariances
% dxi_kf = cell2mat( arrayfun(@(kk) SE2.adjoint( X_kf(:, :, kk)) * ...
%     dxi_kf( :, kk), idx_range, 'UniformOutput', false));
% P_kf    = reshape( cell2mat( arrayfun(@(kk) SE2.adjoint( X_kf(:, :, kk)) * P_kf( :, :, kk)...
%     * SE2.adjoint( X_kf( :, :, kk))', idx_range, 'UniformOutput', false)), 3, 3, []);
% dxi_batch = cell2mat( arrayfun(@(kk) SE2.adjoint( X_batch(:, :, kk)) * ...
%     dxi_batch( :, kk), idx_range, 'UniformOutput', false));
% P_batch = reshape( cell2mat( arrayfun(@(kk) SE2.adjoint( X_batch(:, :, kk)) * P_batch( :, :, kk)...
%     * SE2.adjoint( X_batch( :, :, kk))', idx_range, 'UniformOutput', false)), 3, 3, []);

% Colors
%   Filter
col_kf    = matlabColors( 'orange');
col_batch = matlabColors( 'blue');
%   Error plots
for kk = 1 : 3
    subplot(3, 1, kk);
    hold all; grid on;
    if plt_kf
        %   Filter
        plot( time', dxi_kf(kk,:), 'LineWidth', 1.5, 'DisplayName', 'L-InEKF', ...
            'Color', col_kf);
        plot(  time, - 3 * sqrt( squeeze( P_kf( kk, kk, idx_range))), '-',...
            'LineWidth', 1.5, 'Color', col_kf, 'HandleVisibility', 'off');
        plot(  time, 3 * sqrt( squeeze( P_kf( kk, kk, idx_range))), '-',...
            'LineWidth', 1.5, 'Color', col_kf, 'HandleVisibility', 'off');
    end
    
    if plt_batch
        %   Batch
        plot( time', dxi_batch(kk,:), 'LineWidth', 1.5, 'DisplayName', 'Batch', ...
            'Color', col_batch);
        plot(  time, 3 * sqrt( squeeze( P_batch( kk, kk, idx_range))), '-.',...
            'LineWidth', 1.5, 'Color', col_batch, 'HandleVisibility', 'off');
        plot(  time, - 3 * sqrt( squeeze( P_batch( kk, kk, idx_range))), '-.',...
            'LineWidth', 1.5, 'Color', col_batch, 'HandleVisibility', 'off');
    end
    
    % y-labels
    if kk == 1
        legend({'L-InEKF', 'Batch'}, 'Interpreter', 'latex', 'FontSize', 14);
        ylabel( sprintf('$\\delta \\xi_{%i}$ [rad]', kk),...
                        'Interpreter', 'latex', 'FontSize', 14);
    else
        ylabel( sprintf('$\\delta \\xi_{%i}$ [m]', kk),...
                    'Interpreter', 'latex', 'FontSize', 14);
    end
end
xlabel( '$t_{k}$ [s]', 'Interpreter', 'latex', 'FontSize', 14);

% Plot trajectory
figure; 
plotMlgPose( X_gt_states, '-', matlabColors('grey'));
hold on;
if plt_kf
    plotMlgPose( X_kf_states, '-.', col_kf);
end
if plt_batch
    plotMlgPose( X_batch_states, '-.', col_batch);
end
legend({'Ground truth', 'L-InEKF', 'Batch'}, 'Interpreter', 'latex', 'FontSize', 14);