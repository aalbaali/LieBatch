function [ X_batch, infm_batch] = batchOptimizationSE2(struct_prior, ...
  struct_vel, struct_gyro, struct_gps, struct_lc, t_sim, X_initial, batch_params,...
  optim_params)
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
    %   struct_lc       :   struct
    %       LC measurements, where n_lc is the number of loop closures
    %       Measurement model is of the form
    %         [-]_k = T_1 \ T_2
    %       Contains
    %           meas        :   [ 3 x 3 x n_lc ] SE2 measurements
    %           cov         :   [ 3 x 3 x n_lc ] SPD matrices
    %           idx         :   [ n_lc x 2 ]     LC indices ( [idx_1, idx_2])
    %           time        :   [ n_lc x 2 ]     Time at the LCs
    %   t_sim       :   [ 1 x K] double
    %       Simulation time steps
    %   X_initial   :   [ 3 x 3 x K]
    %       Optimization initial guess.
    %   batch_params  :   struct
    %     Batch parameters (include_lc, include_gps)
    %   optim_params  :   struct
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
    %   Loop closures
    idx_lc = struct_lc.idx;
    t_lc   = struct_lc.time;    
    y_lc   = struct_lc.mean;
    n_lc   = size( idx_lc, 1);
    
    % Compte covariances/weight function
    Sigma = computeCovarianceErrorFunction( X_initial, struct_prior, ...
      struct_vel, struct_gyro, struct_gps, struct_lc, t_sim, batch_params);
    
    % Lambda error function
    func_err = @( X) errorFunction( X, struct_prior, struct_vel, ...
      struct_gyro, struct_gps, struct_lc, t_sim, batch_params);
    
    % Cholesky factor
%     R_sigma = chol( Sigma, 'upper');
    L_sigma = chol( Sigma, 'lower');
    
    % Cost array
    cost_arr = [];
    % Armijo params
    beta = optim_params.armijo_beta;
    c1   = optim_params.armijo_c1;
    
    % States at iteration j
    X_j = X_initial;    
    for lv1 = 1 : optim_params.max_iterations
        % Compute error function
        [ e, J] = func_err( X_j);
        % Compute search direction
        switch lower( optim_params.lin_solver)
            case 'qr'
                % Reorder vars                                
                % Solve using QR decomposition
                [ ~, R] = qr( L_sigma \ [ J, -e], 0);
                % Search direction
                d_k = R( :, 1 : end - 1) \ R( :, end);

            case 'chol'
                % Get weight matrix (inverse of covariance)                                    
                R = chol( sparse( J' * ( Sigma \ J)));
                d_k = -full( R \ ( R' \ (J' * ( Sigma \e))));
            case '\'
                d_k = -(J' * (Sigma \ J)) \ (J' * (Sigma \ e));        
        end
        % Make sure d_k is a column matrix;
        d_k = d_k(:);
        
        % Check if search direction is a descent direction
        if d_k' * J' * (Sigma \ e) >= 0 && norm( d_k) ~= 0
            warning('Not a descent direction!');            
        end
        
        % Reshape search direction 
        d_k_mat = reshape( d_k , 3, []);
        
        % Armijo back-tracking
        obj_val = (1/2) * e' * (Sigma\e); % Objective function value
        grad_val = d_k' * (J' * ( Sigma \ e));
        for lv2 = 0 : optim_params.armijo_max_iteration-1
            alpha_k = beta^lv2;     
            % Increment X_k as to left-invariant error
            X_j_tmp = reshape( cell2mat( arrayfun( @(kk) X_j( :, :, kk) * ...
                se2alg.expMap( - alpha_k * d_k_mat( :, kk)), 1 : K, ...
                'UniformOutput', false)), 3, 3, []);
            % Compute new error function
            [ e, ~] = func_err( X_j_tmp);
            obj_val_tmp = (1/2)* e' * (Sigma \ e);
            if obj_val_tmp <= obj_val + c1 * alpha_k * grad_val
                % Armijo stopping criterion
                % Solution found
                X_j = X_j_tmp;
                break;
            end                        
        end
        if lv2 == optim_params.armijo_max_iteration-1
            warning('Armijo max iteration reached');            
            X_j = X_j_tmp;            
            break;        
        end
        
        % Update error and Jacobian
        [ e, J] = func_err( X_j_tmp);

        cost_arr = [cost_arr; e' * (Sigma \ e)];
        
        % Stopping criterion
        if norm( J' * (Sigma \ e)) * norm(Sigma, 'fro') <= optim_params.tol_ngrad_stop ...
                        || norm( d_k)/numel(d_k) <= optim_params.tol_ngrad_stop
            fprintf('Batch on r_ba_a converged after %i iterations\n', lv1);            
            successful = true;
            break;
        end
    end
    
    % Prompt warning if optimization didn't converge after max.
    % iterations reached.
    if ~successful && ~ ( norm( J' * (Sigma \ e))*norm(Sigma, 'fro') <  ...
            optim_params.tol_ngrad_stop ...
            || norm( d_k)/numel(d_k) <= optim_params.tol_ngrad_stop)
        warning('Optimization did not converge');
    else
        successful = true;
    end

    % Compute information matrix    
    infm_batch = J' * ( Sigma \ J);
    
    X_batch = X_j;
end

function [ cov_err] = computeCovarianceErrorFunction( X_initial, ...
  struct_prior, struct_vel, struct_gyro, struct_gps, struct_lc, t_sim, batch_params)
    % Computes the covariance on the error function    
    
    dt_func_k = @(kk) t_sim( kk + 1) - t_sim( kk);
    
    % Number of poses
    K = length( t_sim);
    % Number of states per pose (constant value)
    n_x = 3;
        
    % GPS params    
    if batch_params.include_gps
      t_gps   = struct_gps.time;
      n_gps   = 2;
      idx_gps = ceil( t_gps / ( t_sim( 2) - t_sim( 1))); % Assuming frequency is constant
      K_gps   = length( t_gps);
    else
      t_gps = [];
      n_gps = 2;
      idx_gps = [];
      K_gps   = 0;
    end
    
    % LC params
    t_lc   = struct_lc.time;
    idx_lc = struct_lc.idx;
    K_lc   = size( idx_lc, 1);
    
    % Total number of random varaibles
    K_rvs = n_x * K + n_gps * K_gps + n_x * K_lc;
    
    % First, compute the Jacobian of the error function w.r.t. the state        
    jac_prior_n = sparse( [], [], [], 3, K_rvs);
    jac_prior_n( 1 : end, 1 : 3) = speye( n_x);
    %   Odometry         
    jac_odom_w_cell = arrayfun( @(kk) dt_func_k(kk) * speye( n_x), 1 : K - 1, 'UniformOutput', false);
    jac_odom_w = sparse( [], [], [], 3 * (K - 1), K_rvs);
    jac_odom_w(1 : end, 3 + (1:n_x * (K -1))) = blkdiag( jac_odom_w_cell{ :});
    %   GPS
    func_jac_gps_n_k = @(kk) X_initial( 1 : 2, 1 : 2, kk)';
    jac_gps_n_cell = arrayfun(@(kk) func_jac_gps_n_k(kk), idx_gps, 'UniformOutput', false);    
    jac_gps_n = sparse( [], [], [], n_gps * K_gps, K_rvs);
    jac_gps_n( 1 : end, n_x * K + (1 : n_gps * K_gps)) = blkdiag( jac_gps_n_cell{ :});
    %   LCs
    jac_lc_n = sparse( [], [], [], n_x * K_lc, K_rvs);
    jac_lc_n( 1 : end, n_x * K + n_gps * K_gps + (1 : n_x * K_lc)) = speye( n_x * K_lc);
    
    % Augment Jacobians
    jac_err = [ jac_prior_n;
                jac_odom_w;
                jac_gps_n;
                jac_lc_n];
   
            
    % Now, the covariances on the random variables
    %   Prior
    cov_prior = struct_prior.cov;
    %   Odometry (i.e., Q)
    %       First, construct the Q 3D matrix (i.e., cov([theta; vel])
    Q_3d = zeros( 3, 3, K - 1);
    Q_3d( 1, 1, :) = struct_gyro.cov(:,:, 1 : K - 1);
    Q_3d( 2 : 3, 2 : 3, :) = struct_vel.cov( :, :, 1 : K - 1);
    cov_odom  = sparse( blkdiag3d( Q_3d));
    %   GPS
    if batch_params.include_gps
      cov_gps   = sparse( blkdiag3d( struct_gps.cov));
    else
      cov_gps   = [];
    end
    %   LC
    cov_lc    = sparse( blkdiag3d( struct_lc.cov));
    
    % Augment covariances (in the block-diagonal sense)
    cov_rvs = blkdiag( cov_prior, cov_odom, cov_gps, cov_lc);    
    
    % Finally, compute the covariances on the error function
    cov_err = (jac_err * cov_rvs * jac_err');
    
    % Ensure symmetry
    cov_err = (1/2) * ( cov_err + cov_err');    
end

function [ err_val, err_jac] = errorFunction( X, struct_prior, struct_vel, ...
  struct_gyro, struct_gps, struct_lc, t_sim, batch_params)
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
    if batch_params.include_gps
      y_gps   = struct_gps.mean;
      t_gps   = struct_gps.time;
      n_gps   = size( y_gps, 1);
      K_gps   = length( t_gps);
      idx_gps = ceil( t_gps / ( t_sim( 2) - t_sim( 1))); % Assuming frequency is constant
    else
      y_gps   = nan( 2, 0);
      t_gps   = [];
      n_gps   = 2;      
      K_gps   = 0;
      idx_gps = [];
    end    
    %   LC
    y_lc    = struct_lc.mean;
    idx_lc  = struct_lc.idx;
    t_lc    = struct_lc.time;
    K_lc    = size( idx_lc, 1);
    
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
            [ C, r] = SE2.decompose( X( :, :, kk));
            err_gps( :, idx_gps_j) = C' * ( y_k - r);            
            % Jacobians of measurement model 
            %   w.r.t. state
            % NOTE: SEEMS TO BE WORKING WITHOUT THE '-' SIGN. I NEED TO FIGURE OUT
            % WHAT'S HAPPENING
            H_k = [ so2alg.wedge(1) * C' * ( y_k - r), eye( 2)];
            jac_gps_x = [ jac_gps_x;
                    kron( sparse( 1, kk, 1, 1, K), H_k)];
            
            % Increment gps index
            idx_gps_j = idx_gps_j + 1;
        end
    end
    
    %   LC
    %       Error array
    err_lc  = nan( n_x, K_lc);
    %       Jacobian of LC error w.r.t. states
    jac_lc_x = sparse( [], [], [], n_x * K_lc, n_x * K);
    for lv1 = 1 : K_lc
      idx_1 = idx_lc( lv1, 1);
      idx_2 = idx_lc( lv1, 2);
      T_1 = X( :, :, idx_1);
      T_2 = X( :, :, idx_2);
      
      y_lc_j = y_lc( :, :, lv1);
      % Error 
      err_lc( :, lv1) = SE2.Log( T_2 \ T_1 * y_lc_j);
      
      % Jacobian of the error function w.r.t. first pose T_1
      jac_lcj_1 = kron( sparse( 1, idx_1, 1, 1, K), ...
        -SE2.computeJRightInv(err_lc( :, lv1)) * SE2.adjoint(SE2.inverse(y_lc_j)));
      % Jacobian of the error function w.r.t. second pose T_2
      jac_lcj_2 = kron( sparse( 1, idx_2, 1, 1, K), SE2.computeJLeftInv( err_lc( :, lv1) ));
      % Augment Jacobians and add to LC Jacobian w.r.t. all states.
      jac_lc_x( (lv1 - 1) * 3 + (1:3), :) = jac_lcj_1 + jac_lcj_2;
    end
    
    % Augment the error arrays and Jacobians
    err_val = [ err_prior; 
                reshape( err_odom, [], 1); 
                reshape( err_gps, [], 1);
                reshape( err_lc, [], 1)];
    err_jac = [ jac_prior_x; 
                jac_odom_x; 
                jac_gps_x;
                jac_lc_x];    
end