function [xhat, P] = UKF(dynamics, x0, ytilde, tmeas, P0, Q, R, alpha, beta, kappa, mu, obsv_lat, LST, R_obsv)

% Storage of outputs
xhat = zeros(length(tmeas), 6); % Contains history of state estiamtes at kth row
P = zeros(6, 6, length(tmeas)); % Contains history of state estimate covariance

% Establishing dimensions
n = size(x0, 1);
q = size(Q, 1);
m = size(R, 1);
n_meas = size(ytilde, 1);
L = n + q + m;

% Initializing variables
% Initial estimate and covariance
xhat_k = x0;
P_k = P0;

% Augmenting to include noise
xhat_k_aug = [xhat_k; zeros(n, 1); zeros(m, 1)];
P_k_aug = blkdiag(P_k, Q, R);

% Scaling Parameters
lambda = alpha^2 * (L + kappa)- L;
gamma = sqrt(L + lambda);

% UKF Weights
n_sigma = 2*L+1; % Total number of sigma points in UKF

% Array to handle the weights for mean and covariance for each sigma point
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
    % Generating sigma_points
    % Using Cholesky decomposition to obtain sqrt(P) using P = S*S'
    S = chol(P_k_aug, "lower");

    % Generating cloud of sigma points
    Chi_k = [xhat_k_aug, xhat_k_aug + gamma*S, xhat_k_aug - gamma*S];

    % Extracting state variables 
    Chi_x_k = Chi_k(1:n, :);

    % Extracting process noise
    Chi_w_k = Chi_k(n+1:n+q, :);

    % Extracting noise variables
    Chi_v_k = Chi_k(n+q+1:end, :);

    % Propogating state of each sigma point from previous step
    % If this is the first step, there is no propogation 
    % Instead, moving directly to update step
    if k == 1
        Chi_x_prop = Chi_x_k;
    else
        Chi_x_prop = zeros(n, n_sigma);
        options = odeset('RelTol', 1E-8, 'AbsTol', 1E-10);
        twobody_ode = @(t, y) dynamics(t, y, mu);
        for l = 1:n_sigma
            [~, simulated_trajectory] = ode45(twobody_ode, tmeas(k-1:k), Chi_x_k(:, l), options);
            Chi_x_prop(:, l) = simulated_trajectory(end, :)' + Chi_w_k(:, l);
        end
    end
    
    % Computing the mean of the sigma points to produce xhat_minus
    % This multiplies the n-th column of Chi_x_prop by the n-th element of
    % W_mean with dimensions 
    % -> (6x2L+1) * (2L+1,1) = (6, 1)
    xhat_minus = Chi_x_prop * W_mean(:);

    % Computing the covariance of the sigma points to produce P_minus
    P_minus = zeros(n, n);
    for l = 1:n_sigma
        P_minus = P_minus + W_cov(l) * (Chi_x_prop(:, l) - xhat_minus) * (Chi_x_prop(:, l) - xhat_minus)';
    end

    % Computing expected observations at each point cloud
    Gamma_k = zeros(m, n_sigma);
    for l = 1:n_sigma
        % Converting from state to measurements
        r_sigma = Chi_x_prop(1:3, l);
        rho_inertial_l = inertial_range(r_sigma', R_obsv, obsv_lat, LST(k));
        rho_obsv_l = observer_range(rho_inertial_l, obsv_lat, LST(k));
        y_sigma = horizontal_coordinates(rho_obsv_l);

        % Adding noise
        Gamma_k(:, l) = y_sigma' + Chi_v_k(:, l);
    end

    % Obtaining yhat_minus as the mean of the expected observations
    yhat_minus = Gamma_k * W_mean(:);

    % Measurement update equations
    % P^eyey is the covaraince of the expected observations above
    % P^exey is the cross-covariance
    P_yy = zeros(m, m);
    P_xy = zeros(n, m);
    for l = 1:n_sigma
        x_err = Chi_x_prop(:, l) - xhat_minus;
        y_err = (Gamma_k(:, l) - yhat_minus);
        
        P_yy = P_yy + W_cov(l) * (y_err * y_err');
        P_xy = P_xy + W_cov(l) * (x_err * y_err');
    end

    % Calculating gain
    K_k = P_xy / P_yy;

    % updating state estimate and covariance
    xhat_k = xhat_minus + K_k * (ytilde(k, :)' - yhat_minus);
    P_k = P_minus - K_k * P_yy * K_k';

    % Storing state estimate and covariance in history
    xhat(k, :) = xhat_k';
    P(:, :, k) = P_k;

    % Updating augmented values for iteration
    xhat_k_aug = [xhat_k; zeros(n, 1); zeros(m, 1)];
    P_k_aug = blkdiag(P_k, Q, R);
end
end