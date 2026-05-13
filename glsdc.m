function [estimate, Lambda] = glsdc(dynamics, guess, measurements, tol, R, tmeas, maxiter, R_obsv, LST, obsv_lat)

% Initializing relevant values
W = inv(R); % 3x3
Phi0 = eye(6);

% Cost for iteration convergence
old_cost = inf;

for n = 1:maxiter

    % Lambda and N for normal equations
    Lambda = zeros(6);
    N = zeros(6, 1);

    % Initial state including r, v, and STM
    initial_state_guess = [guess; Phi0(:)];

    % Integrate dynamics and STM
    options = odeset('RelTol',1E-8, 'AbsTol',1E-10);
    [~, state_estimates] = ode45(dynamics, tmeas, initial_state_guess, options);

    % Converting to inertial, then obsv, then horizontal coordinates to
    % arrived at measurement estimate
    rho = inertial_range(state_estimates(:, 1:3), R_obsv, obsv_lat, LST);
    rho_obsv = observer_range(rho, obsv_lat, LST);
    estimated_measurements = horizontal_coordinates(rho_obsv);

    % Using estimated measurements to calculate error
    err = measurements - estimated_measurements; % 11 x 3

    % Going through and calculating Lambda and N
    new_cost = 0;
    for k = 1:length(tmeas)
        % Calculate H_k = dh/dx * Phi
        % Expressions for dh/dx
        % dh/dx = [dh/dr, dh/dv]; dh/dv = 0
        dh_dx = zeros(3, 6);
        rho_u = rho_obsv(k, 1);
        rho_e = rho_obsv(k, 2);
        rho_n = rho_obsv(k, 3);
        rho_mag = estimated_measurements(k, 1);
        LST_k = LST(k);

        % d(rho) / dr
        drho_dx = (rho_u*cos(obsv_lat)*cos(LST_k) - rho_e*sin(LST_k) - rho_n*sin(obsv_lat)*cos(LST_k))/rho_mag;
        drho_dy = (rho_u*cos(obsv_lat)*sin(LST_k) + rho_e*cos(LST_k) - rho_n*sin(obsv_lat)*sin(LST_k))/rho_mag;
        drho_dz = (rho_u*sin(obsv_lat) + rho_n*cos(obsv_lat))/rho_mag;

        % d(az) / dr
        daz_dx = (rho_e*sin(obsv_lat)*cos(LST_k) - rho_n*sin(LST_k)) / (rho_n^2 + rho_e^2);
        daz_dy = (rho_e*sin(obsv_lat)*sin(LST_k) + rho_n*cos(LST_k)) / (rho_n^2 + rho_e^2);
        daz_dz = -rho_e*cos(obsv_lat) / (rho_n^2 + rho_e^2);

        % d(el) / dr
        del_dx = (rho_mag*cos(obsv_lat)*cos(LST_k) - rho_u*drho_dx) / ( rho_mag*sqrt(rho_mag^2 - rho_u^2) );
        del_dy = (rho_mag*cos(obsv_lat)*sin(LST_k) - rho_u*drho_dy) / ( rho_mag*sqrt(rho_mag^2 - rho_u^2) );
        del_dz = (rho_mag*sin(obsv_lat) - rho_u*drho_dz) / ( rho_mag*sqrt(rho_mag^2 - rho_u^2) );

        % Combining
        dh_dx(1:3, 1:3) = [drho_dx, drho_dy, drho_dz; daz_dx, daz_dy, daz_dz; del_dx, del_dy, del_dz];

        % Extracting Phi
        Phi_k = reshape(state_estimates(k, 7:42), 6, 6);

        % H_k
        H_k = dh_dx * Phi_k;

        % Building Lambda and N
        Lambda = Lambda + H_k' * W * H_k;
        N = N + H_k' * W * err(k, :)';

        % Adding to cost function
        new_cost = new_cost + err(k, :) * W * err(k, :)';
    end

    % Check for convergence
    if n > 1 && abs(new_cost - old_cost)/new_cost < tol
        break;
    end

    % Solve for delta_x
    delta_x = Lambda \ N;

    % Updating guess
    guess = guess + delta_x;
    old_cost = new_cost;
    
end

estimate = guess;
end