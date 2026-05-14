%% Initial Conditions and Constants
clear; clc; close all;
r0 = [7000; 1000; 200]; % km
v0 = [4; 7; 2]; % km/s
y0 = [r0; v0]; % Combined initial state

% Physical
R_obsv = 6371;% Radius of Earth in km
omega_E = 7.2921159E-5; % rad/s
mu = 398600.4415; % Gravitational paramter of earth

% Positional Constants
obsv_lat = deg2rad(5); % Observer latitude
LST0 = deg2rad(10); % Observer local siderial time

% tspan for trajectory integration
tspan = [0, 3000];

% time values for measurements
delta_t = 5;
tmeas = 0:delta_t:3000;

% Computing LST at tmeas
LST = LST0 + omega_E * tmeas;

% Noise Covariance
sigma_rho = 1; % km
sigma_az = deg2rad(0.01); % rad
sigma_el = deg2rad(0.01); % rad
R = diag([sigma_rho^2, sigma_az^2, sigma_el^2]);

% Delimeter for text/section outputs
delim_eq = repelem('=', 70);
delim_dash = repelem('-', 70);

%% Generating Measurements
% Simulating true motion
func = @(t, y) twobody_J2(t, y, mu);
options = odeset('RelTol',1E-10, 'AbsTol',1E-12);
[~, measured_states] = ode45(func, tmeas, y0, options);

% Generating inertial range vector
rho = inertial_range(measured_states(:, 1:3), R_obsv, obsv_lat, LST);

% Range vector in observer frame
rho_obsv = observer_range(rho, obsv_lat, LST);

% Collecting measurements
h_cords = horizontal_coordinates(rho_obsv);

% Adding Noise
measurements = generate_measurements(h_cords, R, length(tmeas));

%% Plotting measurements vs true trajectory in observer frame
[t, true_motion] = ode45(func, tspan, y0, options);
LST_true = LST0 + omega_E*t';
rho_true_motion = inertial_range(true_motion(:, 1:3), R_obsv, obsv_lat, LST_true);
rho_obsv_true_motion = observer_range(rho_true_motion, obsv_lat, LST_true);
h_cords_true_motion = horizontal_coordinates(rho_obsv_true_motion);


%% Batched Least Squares
r0hat = [6990; 1; 1];
v0hat = [1; 1; 1];
x0hat = [r0hat; v0hat];

% BLS and GLSDC use a short arc only -- Gauss-Newton diverges from the cold
% guess over the full 3000 s arc (the linearization is invalid that far out).
N_short     = 30;
tmeas_short = tmeas(1:N_short);
meas_short  = measurements(1:N_short, :);
LST_short   = LST(1:N_short);

dynamics = @(t, y) twobody_STM(t, y, mu);
[bls_estimate, Lambda_bls, final_cost] = bls(dynamics, x0hat, meas_short, R, ...
                                             tmeas_short, R_obsv, LST_short, obsv_lat);
P_bls = inv(Lambda_bls);

%% GLSDC

maxiter = 10;
tol = 1E-3;
dynamics = @(t, y) twobody_STM(t, y, mu);
[glsdc_estimate, Lambda_glsdc] = glsdc(dynamics, x0hat, meas_short, tol, R, ...
                                       tmeas_short, maxiter, R_obsv, LST_short, obsv_lat);
P_glsdc = inv(Lambda_glsdc);

%% RLS w FF

lambda_ff = 0.99;
dynamics = @(t, y) twobody_STM(t, y, mu);
[rls_estimate, Lambda_rls] = rls_ff(dynamics, x0hat, measurements, lambda_ff, R, ...
                                    tmeas, maxiter, R_obsv, LST, obsv_lat);
P_rls = inv(Lambda_rls);

%% Monte Carlo GLSDC Method
nruns = 1000;
mc_results = zeros(nruns, 6);

measurements_orig = measurements;

for n = 1:nruns
    measurements = generate_measurements(h_cords, R, length(tmeas));
    [mc_estimate, Lambda_mc] = glsdc(dynamics, x0hat, measurements(1:N_short, :), tol, R, ...
                                     tmeas_short, maxiter, R_obsv, LST_short, obsv_lat);
    mc_results(n, :) = mc_estimate';
end

% Calculating mean and covariance of estimates from monte carlo runs
monte_carlo_avg = mean(mc_results)';
P_mc = cov(mc_results);

measurements = measurements_orig;

% Error statistics
mc_err = mc_results - y0';
mc_mean = mean(mc_err);
mc_median = median(mc_err);
mc_rmse = sqrt(mean(mc_err.^2));
mc_rmse_pos = sqrt(mean(sum(mc_err(:,1:3).^2, 2)));

%% EKF 

% Initial covariance — square of order of magnitude of initial guess error
P0 = diag([1E6, 1E6, 1E6, 1E2, 1E2, 1E2]);

% Process noise covariance — high confidence in dynamics
% Q = eye(6) * 1E-8;
Q = diag([1e-9, 1e-9, 1e-9, 1e-6, 1e-6, 1e-6]);

% Run EKF
[xhat_ekf, P_ekf] = EKF(x0hat, measurements, tmeas, P0, Q, R, mu, obsv_lat, LST, R_obsv);

% Compute error vs truth (measured_states from ode45 above)
ekf_error = measured_states - xhat_ekf;

% Extract 3-sigma bounds
sigma_bounds = zeros(length(tmeas), 6);
for k = 1:length(tmeas)
    sigma_bounds(k, :) = 3 * sqrt(diag(P_ekf(:, :, k)));
end

% Plot and log
plot_ekf = false;
if plot_ekf
    plotting_ekf;
end

%% UKF

P0 = diag([1E6, 1E6, 1E6, 1E2, 1E2, 1E2]);
% Q  = eye(6) * 1E-8;
Q = diag([1e-9, 1e-9, 1e-9, 1e-6, 1e-6, 1e-6]);

% UKF tuning per Project 5 statement
alpha = 1E-3;
beta  = 2;     % optimal for Gaussian
kappa = 0;

% Run UKF — note @twobody is passed in so the local function can re-wrap it with mu
[xhat_ukf, P_ukf] = UKF(@twobody, x0hat, measurements, tmeas, P0, Q, R, ...
                        alpha, beta, kappa, mu, obsv_lat, LST, R_obsv);

ukf_error = measured_states - xhat_ukf;

sigma_bounds_ukf = zeros(length(tmeas), 6);
for k = 1:length(tmeas)
    sigma_bounds_ukf(k, :) = 3 * sqrt(diag(P_ukf(:, :, k)));
end

plot_ukf = false;
if plot_ukf
    plotting_ukf;
end

%% Warm starting EKF and UKF

% Perform GLSDC for first ~30 time steps
N_warm = 40;
tmeas_warm = tmeas(1:N_warm);
meas_warm  = measurements(1:N_warm, :);
LST_warm   = LST(1:N_warm);

[x0_warm, Lambda_warm] = glsdc(dynamics, x0hat, meas_warm, tol, R, ...
                               tmeas_warm, maxiter, R_obsv, LST_warm, obsv_lat);
P0_warm = inv(Lambda_warm);

prop_opts = odeset('RelTol', 1E-10, 'AbsTol', 1E-12);
[~, x0_warm_prop] = ode45(@(t,y) twobody(t,y,mu), ...
                          [tmeas(1), tmeas(N_warm+1)], x0_warm, prop_opts);
x_warm_start = x0_warm_prop(end, :)';

% Inflate slightly to avoid overconfidence from the batch
P0_warm = diag([1, 1, 1, 1E-2, 1E-2, 1E-2]);

% Run filters from t = tmeas(N_warm) onward
tmeas_post = tmeas(N_warm+1:end);
meas_post  = measurements(N_warm+1:end, :);
LST_post   = LST(N_warm+1:end);

[xhat_ekf_warm, P_ekf_warm] = EKF(x_warm_start, meas_post, tmeas_post, P0_warm, Q, R, ...
                                  mu, obsv_lat, LST_post, R_obsv);

[xhat_ukf_warm, P_ukf_warm] = UKF(@twobody, x_warm_start, meas_post, tmeas_post, P0_warm, Q, R, ...
                                  alpha, beta, kappa, mu, obsv_lat, LST_post, R_obsv);

%% Trajectory RMSE comparison

options = odeset('RelTol', 1E-10, 'AbsTol', 1E-12);

% Propagate the initial-state estimators forward over tmeas
[~, bls_traj]   = ode45(@(t,y) twobody(t,y,mu), tmeas, bls_estimate,   options);
[~, rls_traj]   = ode45(@(t,y) twobody(t,y,mu), tmeas, rls_estimate,   options);
[~, glsdc_traj] = ode45(@(t,y) twobody(t,y,mu), tmeas, glsdc_estimate, options);

% EKF and UKF already give full state histories
% xhat_ekf, xhat_ukf are already N x 6

% Truth: measured_states (already integrated in main.m over tmeas)
[bls_comp,   bls_pos,   bls_vel]   = trajectory_rmse(bls_traj,   measured_states);
[rls_comp,   rls_pos,   rls_vel]   = trajectory_rmse(rls_traj,   measured_states);
[glsdc_comp, glsdc_pos, glsdc_vel] = trajectory_rmse(glsdc_traj, measured_states);
[ekf_comp,   ekf_pos,   ekf_vel]   = trajectory_rmse(xhat_ekf,   measured_states);
[ukf_comp,   ukf_pos,   ukf_vel]   = trajectory_rmse(xhat_ukf,   measured_states);

% Assemble into a table
estimator   = {'BLS'; 'RLS-FF'; 'GLSDC'; 'EKF'; 'UKF'};
rmse_pos    = [bls_pos;   rls_pos;   glsdc_pos;   ekf_pos;   ukf_pos];
rmse_vel    = [bls_vel;   rls_vel;   glsdc_vel;   ekf_vel;   ukf_vel];
rmse_table  = table(estimator, rmse_pos, rmse_vel);
disp(rmse_table)

%% Warm-started filter RMSE (GLSDC-seeded)
% Warm filters cover only tmeas(N_warm+1:end), so slice truth to match.
measured_states_post = measured_states(N_warm+1:end, :);

[ekf_warm_comp, ekf_warm_pos, ekf_warm_vel] = trajectory_rmse(xhat_ekf_warm, measured_states_post);
[ukf_warm_comp, ukf_warm_pos, ukf_warm_vel] = trajectory_rmse(xhat_ukf_warm, measured_states_post);

% Cold EKF/UKF scored over the SAME tail window, for an apples-to-apples
% warm-vs-cold comparison (removes the window-length confound: the warm
% filters skip the early transient by construction).
xhat_ekf_tail = xhat_ekf(N_warm+1:end, :);
xhat_ukf_tail = xhat_ukf(N_warm+1:end, :);
[ekf_tail_comp, ekf_tail_pos, ekf_tail_vel] = trajectory_rmse(xhat_ekf_tail, measured_states_post);
[ukf_tail_comp, ukf_tail_pos, ukf_tail_vel] = trajectory_rmse(xhat_ukf_tail, measured_states_post);

% Comparison table over the tail window
estimator_tail = {'EKF (cold, tail)'; 'EKF (warm)'; 'UKF (cold, tail)'; 'UKF (warm)'};
rmse_pos_tail  = [ekf_tail_pos; ekf_warm_pos; ukf_tail_pos; ukf_warm_pos];
rmse_vel_tail  = [ekf_tail_vel; ekf_warm_vel; ukf_tail_vel; ukf_warm_vel];
rmse_table_warm = table(estimator_tail, rmse_pos_tail, rmse_vel_tail);
disp(rmse_table_warm)