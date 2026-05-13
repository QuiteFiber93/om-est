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
delta_t = 10;
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
options = odeset('RelTol',1E-10, 'AbsTol',1E-12);
[~, measured_states] = ode89(func, tmeas, y0, options);

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

% Figure Plotting true trajectory
% Creating subplots for x, y, and z
fig = tiledlayout(3, 1);
ax1 = nexttile;
plot(ax1, t, true_motion(:, 1))
ylabel('x(t) (km)')

ax2 = nexttile;
plot(ax2, t, true_motion(:, 2))
ylabel('y(t) (km)')

ax3 = nexttile;
plot(ax3, t, true_motion(:, 3))
ylabel('z(t) (km)')
xlabel('t (s)')
linkaxes([ax1, ax2, ax3], 'x')
title(fig, 'Earth Fixed Position');

exportgraphics(gcf, 'Images/body_fixed_pos.png','Resolution',300)

% Plotting x, y, and z velocities in subplots and saving the figure
fig = tiledlayout(3, 1);
ax1 = nexttile;
plot(ax1, t, true_motion(:, 4))
ylabel('xdot(t) (km/s)')

ax2 = nexttile;
plot(ax2, t, true_motion(:, 5))
ylabel('ydot(t) (km/s)')

ax3 = nexttile;
plot(ax3, t, true_motion(:, 6))
ylabel('zdot(t) (km/s)')
xlabel('t (s)')
title(fig, 'Earth Fixed Velocity');

exportgraphics(gcf, 'Images/body_fixed_vel.png','Resolution',300)

figure
% converting true trajectory from spherical (horizontal) frame to cartesian
rho_true = h_cords_true_motion(:, 1);
az_true = h_cords_true_motion(:, 2);
el_true = h_cords_true_motion(:, 3);
[x_obsv, y_obsv, z_obsv] = sph2cart(az_true, el_true, rho_true);
plot3(x_obsv, y_obsv, z_obsv, 'DisplayName','True Trajectory');
hold on
[x_meas, y_meas, z_meas] = sph2cart(measurements(:, 2), measurements(:, 3), measurements(:, 1));
scatter3(x_meas, y_meas, z_meas, 10, 'filled', 'DisplayName', 'Measured Positions');
hold off
title('True Trajectory vs Measured States')
grid()
legend('Location','northeast')
xlabel('x (km)')
ylabel('y (km)')
zlabel('z (km)')

% Saving plot
exportgraphics(gca, 'Images/measurements.png','Resolution',300)

% Plotting x, y, z trajectory relative to the observer vs time and
% overlaying measurements
figure;
ax1 = subplot(3, 1, 1);
hold on
plot(ax1, t, x_obsv, 'DisplayName','True')
scatter(ax1, tmeas, x_meas, 10, 'filled', 'DisplayName','Measured')
hold off
ylabel('x(t) (km)')
legend('Location','southeast')

ax2 = subplot(3, 1, 2);
hold on
plot(ax2, t, y_obsv, 'DisplayName','True')
scatter(ax2, tmeas, y_meas, 10, 'filled', 'DisplayName','Measured')
hold off
ylabel('y(t) (km)')
legend('Location','southeast')

ax3 = subplot(3, 1, 3);
hold on
plot(ax3, t, z_obsv, 'DisplayName','True')
scatter(ax3, tmeas, z_meas, 10, 'filled', 'DisplayName','Measured')
hold off
ylabel('z(t) (km)')
xlabel('t (s)')
legend('Location','southeast')

sgtitle('Position of Satellite Relative to Observer')
exportgraphics(gcf, 'Images/cartesian_pos_measurements.png', 'Resolution',300)

% Creating a figure which shows the altitude and azimuth of the satellite
% overhead from frame of observer
figure
polarplot(az_true, 90 - rad2deg(el_true), 'DisplayName','True Trajectory')
hold on
polarscatter(measurements(:, 2), 90 - rad2deg(measurements(:, 3)), 7, 'filled', 'DisplayName','Measurements')
hold off
title('Polar Plot for Altitude and Azimuth')
rlim([0, 90])
ax = gca;
ax.RTick = 0:30:90;
ax.RTickLabel = {'90°', '60°', '30°', '0°'};
exportgraphics(gcf, 'Images/alt_az_plot.png', 'Resolution',300)

% Plotting Range vs Range Measurements
figure
plot(t, rho_true, 'DisplayName','True Range')
hold on
scatter(tmeas, measurements(:, 1), 10, 'filled', 'DisplayName','Range Measurements')
hold off
title('True Range vs Range Measurements')
xlabel('t (s)')
ylabel('Range (km)')
legend()
exportgraphics(gcf, 'Images/range_plot.png', 'Resolution',300)

%% Batched Least Squares
r0hat = [6990; 1; 1];
v0hat = [1; 10; 1];
x0hat = [r0hat; v0hat];

dynamics = @(t, y) twobody_STM(t, y, mu);
[bls_estimate, Lambda] = bls(dynamics, x0hat, measurements, R, tmeas, R_obsv, LST, obsv_lat);
P_bls = inv(Lambda);

%% GLSDC

maxiter = 10;
tol = 1E-3;
dynamics = @(t, y) twobody_STM(t, y, mu);
[glsdc_estimate, Lambda] = glsdc(dynamics, x0hat, measurements, tol, R, tmeas, maxiter, R_obsv, LST, obsv_lat);
P_glsdc = inv(Lambda);

%% Monte Carlo GLSDC Method
nruns = 1000;
mc_results = zeros(nruns, 6);

for n = 1:nruns
    measurements = generate_measurements(h_cords, R, length(tmeas));
    [mc_estimate, Lambda] = glsdc(dynamics, x0hat, measurements, tol, R, tmeas, maxiter, R_obsv, LST, obsv_lat);
    mc_results(n, :) = mc_estimate';
end

% Calculating mean and covariance of estimates from monte carlo runs
monte_carlo_avg = mean(mc_results)';
P_mc = cov(mc_results);

%% EKF 

% Initial covariance — square of order of magnitude of initial guess error
P0 = diag([1E6, 1E6, 1E6, 1E2, 1E2, 1E2]);

% Process noise covariance — high confidence in dynamics
Q = eye(6) * 1E-8;

% Run EKF
[xhat_ekf, P_ekf] = EKF(x0hat, measurements, tmeas, P0, Q, R, mu, obsv_lat, LST, R_obsv);

% Compute error vs truth (measured_states from ode89 above)
ekf_error = measured_states - xhat_ekf;

% Extract 3-sigma bounds
sigma_bounds = zeros(length(tmeas), 6);
for k = 1:length(tmeas)
    sigma_bounds(k, :) = 3 * sqrt(diag(P_ekf(:, :, k)));
end

% Plot and log
plot_ekf = true;
if plot_ekf
    plotting_ekf;
end

%% UKF

P0 = diag([1E6, 1E6, 1E6, 1E2, 1E2, 1E2]);
Q  = eye(6) * 1E-8;

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

plot_ukf = true;
if plot_ukf
    plotting_ukf;
end

%% Warm starting EKF and UKF

% Perform GLSDC for first ~30 time steps
N_warm = 30;
tmeas_warm = tmeas(1:N_warm);
meas_warm  = measurements(1:N_warm, :);
LST_warm   = LST(1:N_warm);

[x0_warm, Lambda_warm] = glsdc(dynamics, x0hat, meas_warm, tol, R, ...
                               tmeas_warm, maxiter, R_obsv, LST_warm, obsv_lat);
P0_warm = inv(Lambda_warm);

% Inflate slightly to avoid overconfidence from the batch
P0_warm = P0_warm + diag([1, 1, 1, 1E-2, 1E-2, 1E-2]);

% Run filters from t = tmeas(N_warm) onward
tmeas_post = tmeas(N_warm:end);
meas_post  = measurements(N_warm:end, :);
LST_post   = LST(N_warm:end);

[xhat_ekf_warm, P_ekf_warm] = EKF(x0_warm, meas_post, tmeas_post, P0_warm, Q, R, ...
                                  mu, obsv_lat, LST_post, R_obsv);

[xhat_ukf_warm, P_ukf_warm] = UKF(@twobody, x0_warm, meas_post, tmeas_post, P0_warm, Q, R, ...
                                  alpha, beta, kappa, mu, obsv_lat, LST_post, R_obsv);