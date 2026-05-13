function [xhat, P] = EKF(x0, ytilde, tmeas, P0, Q, R, mu, obsv_lat, LST, R_obsv)

n_meas = length(tmeas);

% Creating variables to store results
xhat = zeros(n_meas, 6);
P = zeros(6, 6, n_meas);

% Initializing
xhat_minus = x0;
P_minus = P0;

% Iterating through measurements
for k = 1:n_meas
    % Computing H_k and expected measurements
    H_k = compute_H(xhat_minus, obsv_lat, LST(k), R_obsv);

    r_inertial = xhat_minus(1:3);
    rho_inertial = inertial_range(r_inertial', R_obsv, obsv_lat, LST(k));
    rho_obsv = observer_range(rho_inertial, obsv_lat, LST(k));
    y_prediction = horizontal_coordinates(rho_obsv)';

    % Computing Kalman Gain
    K_k = P_minus * H_k' / (H_k * P_minus * H_k' + R);

    % Updating state and covariance estimates
    xhat_plus = xhat_minus + K_k * (ytilde(k, :)' - y_prediction);
    P_plus = (eye(6) - K_k*H_k) * P_minus * (eye(6) - K_k*H_k)' + K_k * R * K_k';

    % Storing results
    xhat(k, :) = xhat_plus';
    P(:, :, k) = P_plus;

    % Propagating to the next time step
    if k < n_meas
        % Defining time span and ODE options
        tspan = [tmeas(k), tmeas(k+1)];
        options = odeset('RelTol',1E-8, 'AbsTol',1E-10);

        % Propagating state and covariance
        initial_state = [xhat_plus; P_plus(:)];
        [~, state_prop] = ode45(@(t, y) covar_dynamics(t, y, Q, mu), tspan, initial_state, options);
        next_step = state_prop(end, :)';
        xhat_minus = next_step(1:6);
        P_minus = reshape(next_step(7:42), 6, 6);

    end
end
end
