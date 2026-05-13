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
delta_t = 100;
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
func = @(t, y) twobody(t, y, mu);
options = odeset('RelTol',1E-8, 'AbsTol',1E-10);
[~, measured_states] = ode45(func, tmeas, y0, options);

% Generating inertial range vector
rho = inertial_range(measured_states(:, 1:3), R_obsv, obsv_lat, LST);

% Range vector in observer frame
rho_obsv = observer_range(rho, obsv_lat, LST);
%
% % Collecting measurements
h_cords = horizontal_coordinates(rho_obsv);
% h_cords = measurement_model(measured_states(:, 1:3), obsv_lat, LST, R_obsv);
% Adding Noise
measurements = generate_measurements(h_cords, R, length(tmeas));

%% Plotting measurements vs true trajectory in observer frame

[t, true_motion] = ode45(func, tspan, y0, options);
LST_true = LST0 + omega_E*t';
rho_true_motion = inertial_range(true_motion(:, 1:3), R_obsv, obsv_lat, LST_true);
rho_obsv_true_motion = observer_range(rho_true_motion, obsv_lat, LST_true);
h_cords_true_motion = horizontal_coordinates(rho_obsv_true_motion);

% converting true trajectory from spherical (horizontal) frame to cartesian
rho_true = h_cords_true_motion(:, 1);
az_true = h_cords_true_motion(:, 2);
el_true = h_cords_true_motion(:, 3);
[x_obsv, y_obsv, z_obsv] = sph2cart(az_true, el_true, rho_true);

% Decides to actually plot
plot_measurements = true;
if plot_measurements
    plotting_measurements;
end

%% Implementing Extended Kalman Filter

% Initial guess
x0_ekf = [6990; 1; 1; 1; 1; 1];
P0 = diag([1E6, 1E6, 1E6, 1E2, 1E2, 1E2]);

% Process noise covariance
Q = eye(6) * 1E-8;

% Run EKF
[xhat_ekf, P_ekf] = EKF(x0_ekf, measurements, tmeas, P0, Q, R, mu, obsv_lat, LST, R_obsv);

% Calculating error
ekf_error = measured_states - xhat_ekf;

% Extracting 3 sigma error bounds
sigma_bounds = zeros(length(tmeas), 6);
for k=1:length(tmeas)
    sigma_bounds(k, :) = 3 * sqrt(diag(P_ekf(:, :, k)));
end

%% Plotting Error For EKF

plot_ekf = true;
if plot_ekf
    plotting_ekf;
    log_filter_results(tmeas, ekf_error, sigma_bounds, delta_t, 'EKF');
end

%% Implementing Unscented Kalman Filter
% Add this section to your main.m after the EKF section

% Initial guess (same as EKF)
x0_ukf = [6990; 1; 1; 1; 1; 1];

% Initial covariance (same as EKF)
P0 = diag([1E6, 1E6, 1E6, 1E2, 1E2, 1E2]);

% Process noise covariance (same as EKF)
Q = eye(6) * 1E-8;

% UKF tuning parameters (as specified in project)
alpha = 1E-3;
beta = 2;       % Optimal for Gaussian distributions
kappa = 0;

% Run UKF
[xhat_ukf, P_ukf] = UKF(@twobody, x0_ukf, measurements, tmeas, P0, Q, R, alpha, beta, kappa, mu, obsv_lat, LST, R_obsv);

% Calculating error
ukf_error = measured_states - xhat_ukf;

% Extracting 3 sigma error bounds
sigma_bounds_ukf = zeros(length(tmeas), 6);
for k = 1:length(tmeas)
    sigma_bounds_ukf(k, :) = 3 * sqrt(diag(P_ukf(:, :, k)));
end


%% Plotting Error For UKF
plot_ukf = true;
if plot_ukf
    plotting_ukf;
    log_filter_results(tmeas, ukf_error, sigma_bounds, delta_t, 'UKF');
end
