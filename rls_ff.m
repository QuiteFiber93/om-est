function [estimate, Lambda] = rls_ff(dynamics, guess, measurements, lambda_ff, R, tmeas, R_obsv, LST, obsv_lat)
% Sequential recursive least-squares with forgetting factor for
% initial-state estimation.
%
% Single forward pass through the measurements. Each measurement is
% mapped back to the initial state via the state transition matrix
% Phi(t_k, t_0), and the information matrix / vector are updated
% recursively with exponential discounting of past information.
%
% Inputs:
%   dynamics     - ODE handle for [x_dot; Phi_dot(:)] (state + STM), 42x1
%   guess        - 6x1 initial-state guess (used to linearize the trajectory)
%   measurements - N x 3 array of (range, az, el) measurements
%   lambda_ff    - forgetting factor in (0, 1]; effective memory ~ 1/(1-lambda)
%   R            - 3x3 measurement noise covariance
%   tmeas        - N x 1 vector of measurement times (tmeas(1) = t_0)
%   R_obsv       - observer position magnitude (Earth radius)
%   LST          - N x 1 vector of local sidereal time at each measurement
%   obsv_lat     - observer geodetic latitude (rad)
%
% Outputs:
%   estimate     - 6x1 final estimate of the initial state
%   Lambda       - 6x6 final information matrix at t_0 (inverse = covariance)

    W = inv(R);
    n_state = 6;

    % Propagate the reference trajectory + STM once, about the initial guess.
    % Sequential RLS does NOT re-linearize; it processes measurements in
    % one forward pass.
    Phi0 = eye(n_state);
    initial_state_guess = [guess; Phi0(:)];
    options = odeset('RelTol', 1E-8, 'AbsTol', 1E-10);
    [~, state_estimates] = ode45(dynamics, tmeas, initial_state_guess, options);

    % Predicted measurements along the reference trajectory.
    rho      = inertial_range(state_estimates(:, 1:3), R_obsv, obsv_lat, LST);
    rho_obsv = observer_range(rho, obsv_lat, LST);
    estimated_measurements = horizontal_coordinates(rho_obsv);
    err = measurements - estimated_measurements;

    % Initialize information matrix and vector at t_0.
    % Lambda_0 = 0 means a non-informative prior on the initial state.
    % (If you have a prior covariance P0, set Lambda = inv(P0) and
    % N = Lambda * (prior_mean - guess) instead.)
    Lambda = zeros(n_state);
    N_vec  = zeros(n_state, 1);

    % Single forward pass: recursive update with forgetting factor.
    for k = 1:length(tmeas)
        % Measurement Jacobian dh/dx at time t_k
        rho_u   = rho_obsv(k, 1);
        rho_e   = rho_obsv(k, 2);
        rho_n   = rho_obsv(k, 3);
        rho_mag = estimated_measurements(k, 1);
        LST_k   = LST(k);

        drho_dx = (rho_u*cos(obsv_lat)*cos(LST_k) - rho_e*sin(LST_k) - rho_n*sin(obsv_lat)*cos(LST_k)) / rho_mag;
        drho_dy = (rho_u*cos(obsv_lat)*sin(LST_k) + rho_e*cos(LST_k) - rho_n*sin(obsv_lat)*sin(LST_k)) / rho_mag;
        drho_dz = (rho_u*sin(obsv_lat) + rho_n*cos(obsv_lat)) / rho_mag;

        daz_dx = (rho_e*sin(obsv_lat)*cos(LST_k) - rho_n*sin(LST_k)) / (rho_n^2 + rho_e^2);
        daz_dy = (rho_e*sin(obsv_lat)*sin(LST_k) + rho_n*cos(LST_k)) / (rho_n^2 + rho_e^2);
        daz_dz = -rho_e*cos(obsv_lat) / (rho_n^2 + rho_e^2);

        del_dx = (rho_mag*cos(obsv_lat)*cos(LST_k) - rho_u*drho_dx) / (rho_mag*sqrt(rho_mag^2 - rho_u^2));
        del_dy = (rho_mag*cos(obsv_lat)*sin(LST_k) - rho_u*drho_dy) / (rho_mag*sqrt(rho_mag^2 - rho_u^2));
        del_dz = (rho_mag*sin(obsv_lat) - rho_u*drho_dz) / (rho_mag*sqrt(rho_mag^2 - rho_u^2));

        dh_dx = zeros(3, n_state);
        dh_dx(1:3, 1:3) = [drho_dx, drho_dy, drho_dz; ...
                           daz_dx,  daz_dy,  daz_dz;  ...
                           del_dx,  del_dy,  del_dz];

        % Map to initial-state sensitivity via the STM
        Phi_k = reshape(state_estimates(k, 7:42), n_state, n_state);
        H_k   = dh_dx * Phi_k;

        % Recursive information update with forgetting factor:
        %   Lambda_k = lambda * Lambda_{k-1} + H_k' * W * H_k
        %   N_k     = lambda * N_{k-1}     + H_k' * W * err_k
        % The forgetting factor exponentially down-weights older
        % information as new measurements come in, giving an effective
        % sliding window of ~1/(1-lambda) measurements.
        Lambda = lambda_ff * Lambda + H_k' * W * H_k;
        N_vec  = lambda_ff * N_vec  + H_k' * W * err(k, :)';

        % Symmetrize for numerical stability
        Lambda = 0.5 * (Lambda + Lambda');
    end

    % Solve the normal equations once at the end for the initial-state
    % correction (the STM has already mapped every measurement back to t_0).
    delta_x  = Lambda \ N_vec;
    estimate = guess + delta_x;
end