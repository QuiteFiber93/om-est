function [estimate, Lambda] = rls_ff(dynamics, guess, measurements, lambda_ff, R, tmeas, maxiter, R_obsv, LST, obsv_lat)
% Recursive least-squares with forgetting factor for initial-state estimation.

    W = inv(R);
    Phi0 = eye(6);
    options = odeset('RelTol', 1E-8, 'AbsTol', 1E-10);

    old_cost = inf;
    tol = 1E-3;   % outer-loop convergence tolerance (matches glsdc.m)

    for n = 1:maxiter

        % Integrate trajectory + STM about the current guess
        initial_state_guess = [guess; Phi0(:)];
        [~, state_estimates] = ode45(dynamics, tmeas, initial_state_guess, options);

        % Predicted measurements along this trajectory
        rho       = inertial_range(state_estimates(:, 1:3), R_obsv, obsv_lat, LST);
        rho_obsv  = observer_range(rho, obsv_lat, LST);
        estimated_measurements = horizontal_coordinates(rho_obsv);

        err = measurements - estimated_measurements;

        % --- Recursive (per-measurement) information accumulation ---
        Lambda = zeros(6);
        N_vec  = zeros(6, 1);
        new_cost = 0;

        for k = 1:length(tmeas)
            % Measurement Jacobian block dh/dx (same form as bls.m / glsdc.m)
            dh_dx   = zeros(3, 6);
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

            dh_dx(1:3, 1:3) = [drho_dx, drho_dy, drho_dz; ...
                               daz_dx,  daz_dy,  daz_dz;  ...
                               del_dx,  del_dy,  del_dz];

            % Map to initial-state sensitivity via the STM
            Phi_k = reshape(state_estimates(k, 7:42), 6, 6);
            H_k   = dh_dx * Phi_k;

            % Recursive update with forgetting factor: discount prior
            % information, then add this measurement's contribution.
            Lambda = 0.5 * (Lambda + Lambda');
            Lambda = lambda_ff * Lambda + H_k' * W * H_k;
            N_vec  = lambda_ff * N_vec  + H_k' * W * err(k, :)';

            new_cost = new_cost + err(k, :) * W * err(k, :)';
        end

        % Outer Gauss-Newton step: solve normal equations, update guess
        delta_x = Lambda \ N_vec;
        guess   = guess + delta_x;

        % Convergence check on the relative cost change (as in glsdc.m)
        if n > 1 && abs(new_cost - old_cost) / new_cost < tol
            break;
        end
        old_cost = new_cost;

    end

    estimate = guess;
end