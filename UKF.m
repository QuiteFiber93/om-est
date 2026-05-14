function [xhat, P] = UKF(dynamics, x0, ytilde, tmeas, P0, Q, R, alpha, beta, kappa, mu, obsv_lat, LST, R_obsv)

% Storage of outputs
xhat = zeros(length(tmeas), 6); % Contains history of state estimates at kth row
P = zeros(6, 6, length(tmeas)); % Contains history of state estimate covariance

% Establishing dimensions
n = size(x0, 1);
m = size(R, 1);
n_meas = size(ytilde, 1);

% Augmented dimension: state + measurement noise only.
% Process noise is additive (consistent with the EKF's Pdot = F P + P F' + Q),
% so it is NOT carried in the augmented state. Instead it enters as an
% explicit Q_d added to the predicted covariance below. The measurement
% noise block IS augmented, per the professor's Fig. 3 formulation.
L = n + m;

% Initializing variables
xhat_k = x0;
P_k = P0;

% Augmenting to include measurement noise
xhat_k_aug = [xhat_k; zeros(m, 1)];
P_k_aug = blkdiag(P_k, R);

% Scaling Parameters
lambda = alpha^2 * (L + kappa) - L;
gamma = sqrt(L + lambda);

% UKF Weights
n_sigma = 2*L + 1; % Total number of sigma points

W_mean = zeros(n_sigma, 1);
W_cov = zeros(n_sigma, 1);

% Weight of 0th sigma point
W_mean(1) = lambda / (L + lambda);
W_cov(1) = lambda / (L + lambda) + (1 - alpha^2 + beta);

% Weight of all other sigma points
for k = 2:n_sigma
    W_mean(k) = 1 / (2 * (L + lambda));
    W_cov(k) = 1 / (2 * (L + lambda));
end

% Loop for sequential estimation
for k = 1:n_meas
    % Generating sigma points
    % Cholesky decomposition for sqrt(P), symmetrized for stability
    P_k_aug = 0.5 * (P_k_aug + P_k_aug');
    S = chol(P_k_aug, "lower");

    % Cloud of sigma points
    Chi_k = [xhat_k_aug, xhat_k_aug + gamma*S, xhat_k_aug - gamma*S];

    % Extracting state variables
    Chi_x_k = Chi_k(1:n, :);

    % Extracting measurement-noise variables
    Chi_v_k = Chi_k(n+1:end, :);

    % Propagating state of each sigma point from previous step.
    % First step has no propagation -- go straight to the update.
    if k == 1
        Chi_x_prop = Chi_x_k;
        Q_d = zeros(n, n);   % no propagation interval, no process noise yet
    else
        % First-order discretization of the continuous-time process noise
        % spectral density Q over the propagation interval. This makes the
        % UKF's per-step process noise match what the EKF accumulates by
        % integrating Pdot = F P + P F' + Q. (Van Loan's method would give
        % the exact Q_d if higher fidelity were needed.)
        dt = tmeas(k) - tmeas(k-1);
        Q_d = Q * dt;

        Chi_x_prop = zeros(n, n_sigma);
        options = odeset('RelTol', 1E-8, 'AbsTol', 1E-10);
        twobody_ode = @(t, y) dynamics(t, y, mu);
        for l = 1:n_sigma
            [~, simulated_trajectory] = ode45(twobody_ode, tmeas(k-1:k), Chi_x_k(:, l), options);
            Chi_x_prop(:, l) = simulated_trajectory(end, :)';
        end
    end

    % Mean of the propagated sigma points -> xhat_minus
    % (6 x n_sigma) * (n_sigma x 1) = (6 x 1)
    xhat_minus = Chi_x_prop * W_mean(:);

    % Covariance of the propagated sigma points, plus additive process noise Q_d
    P_minus = zeros(n, n);
    for l = 1:n_sigma
        P_minus = P_minus + W_cov(l) * (Chi_x_prop(:, l) - xhat_minus) * (Chi_x_prop(:, l) - xhat_minus)';
    end
    P_minus = P_minus + Q_d;

    % Expected observation at each sigma point
    Gamma_k = zeros(m, n_sigma);
    for l = 1:n_sigma
        r_sigma = Chi_x_prop(1:3, l);
        rho_inertial_l = inertial_range(r_sigma', R_obsv, obsv_lat, LST(k));
        rho_obsv_l = observer_range(rho_inertial_l, obsv_lat, LST(k));
        y_sigma = horizontal_coordinates(rho_obsv_l);

        % Measurement noise enters through the augmented v-block
        Gamma_k(:, l) = y_sigma' + Chi_v_k(:, l);
    end

    % Mean of the expected observations
    yhat_minus = Gamma_k * W_mean(:);

    % Measurement-update covariances
    % P_yy: covariance of expected observations
    % P_xy: state/observation cross-covariance
    P_yy = zeros(m, m);
    P_xy = zeros(n, m);
    for l = 1:n_sigma
        x_err = Chi_x_prop(:, l) - xhat_minus;
        y_err = Gamma_k(:, l) - yhat_minus;

        P_yy = P_yy + W_cov(l) * (y_err * y_err');
        P_xy = P_xy + W_cov(l) * (x_err * y_err');
    end

    % Gain
    K_k = P_xy / P_yy;

    % State and covariance update
    xhat_k = xhat_minus + K_k * (ytilde(k, :)' - yhat_minus);
    P_k = P_minus - K_k * P_yy * K_k';

    % Symmetrize the updated covariance for stability
    P_k = 0.5 * (P_k + P_k');

    % Store
    xhat(k, :) = xhat_k';
    P(:, :, k) = P_k;

    % Rebuild augmented quantities for the next iteration
    xhat_k_aug = [xhat_k; zeros(m, 1)];
    P_k_aug = blkdiag(P_k, R);
end
end