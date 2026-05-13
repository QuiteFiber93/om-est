function [estimate, Lambda, final_cost] = bls(dynamics, guess, measurements, R, tmeas, R_obsv, LST, obsv_lat)
% Batch least-squares: GLSDC with a single Gauss-Newton step.
% Linearizes the nonlinear measurement model once around `guess`, solves
% the normal equations, and returns guess + delta_x.

    W = inv(R);
    Phi0 = eye(6);

    Lambda = zeros(6);
    N = zeros(6, 1);

    initial_state_guess = [guess; Phi0(:)];
    options = odeset('RelTol', 1E-8, 'AbsTol', 1E-10);
    [~, state_estimates] = ode45(dynamics, tmeas, initial_state_guess, options);

    rho       = inertial_range(state_estimates(:, 1:3), R_obsv, obsv_lat, LST);
    rho_obsv  = observer_range(rho, obsv_lat, LST);
    estimated_measurements = horizontal_coordinates(rho_obsv);

    err = measurements - estimated_measurements;

    final_cost = 0;
    for k = 1:length(tmeas)
        dh_dx = zeros(3, 6);
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

        dh_dx(1:3, 1:3) = [drho_dx, drho_dy, drho_dz; daz_dx, daz_dy, daz_dz; del_dx, del_dy, del_dz];
        Phi_k = reshape(state_estimates(k, 7:42), 6, 6);
        H_k = dh_dx * Phi_k;

        Lambda     = Lambda     + H_k' * W * H_k;
        N          = N          + H_k' * W * err(k, :)';
        final_cost = final_cost + err(k, :) * W * err(k, :)';
    end

    delta_x = Lambda \ N;
    estimate = guess + delta_x;
end