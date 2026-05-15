%% ===================================================================
%  PORTFOLIO FIGURES
%  Assumes all estimators have been run and the following exist:
%    tmeas, measured_states, measurements, h_cords_true_motion,
%    t, delta_t,
%    ekf_error, ukf_error, sigma_bounds, sigma_bounds_ukf,
%    bls_traj, rls_traj, glsdc_traj, xhat_ekf, xhat_ukf,
%    bls_pos, bls_vel, rls_pos, rls_vel, glsdc_pos, glsdc_vel,
%    ekf_pos, ekf_vel, ukf_pos, ukf_vel,
%    tmeas_post, measured_states_post,
%    xhat_ekf_warm, xhat_ukf_warm,
%    ekf_warm_pos, ekf_warm_vel, ukf_warm_pos, ukf_warm_vel,
%    ekf_tail_pos, ekf_tail_vel, ukf_tail_pos, ukf_tail_vel
%  ===================================================================

%% Figure 1: Measurement Input Data
figure('Position', [100 100 700 550])

subplot(3,1,1)
plot(t, h_cords_true_motion(:,1), 'b', 'LineWidth', 1)
hold on
scatter(tmeas, measurements(:,1), 4, 'r', 'filled', 'MarkerFaceAlpha', 0.4)
hold off
ylabel('Range (km)')
title('Ground Station Measurements vs Truth')
legend('Truth', 'Noisy Measurements', 'Location', 'best')
grid on

subplot(3,1,2)
plot(t, rad2deg(h_cords_true_motion(:,2)), 'b', 'LineWidth', 1)
hold on
scatter(tmeas, rad2deg(measurements(:,2)), 4, 'r', 'filled', 'MarkerFaceAlpha', 0.4)
hold off
ylabel('Azimuth (deg)')
grid on

subplot(3,1,3)
plot(t, rad2deg(h_cords_true_motion(:,3)), 'b', 'LineWidth', 1)
hold on
scatter(tmeas, rad2deg(measurements(:,3)), 4, 'r', 'filled', 'MarkerFaceAlpha', 0.4)
hold off
ylabel('Elevation (deg)')
xlabel('Time (s)')
grid on

exportgraphics(gcf, 'Images/fig1_measurements.png', 'Resolution', 300)


%% Figure 2: EKF vs UKF Position Error with 3-sigma bounds
figure('Position', [100 100 900 600])
labels = {'x', 'y', 'z'};

for i = 1:3
    % EKF (left column)
    subplot(3, 2, 2*i - 1)
    hold on
    plot(tmeas, ekf_error(:, i), 'b', 'LineWidth', 0.8)
    plot(tmeas,  sigma_bounds(:, i), 'r--', 'LineWidth', 1)
    plot(tmeas, -sigma_bounds(:, i), 'r--', 'LineWidth', 1)
    hold off
    ylabel(sprintf('\\Delta%s (km)', labels{i}))
    grid on
    if i == 1
        title('EKF')
        legend('Error', '3\sigma', 'Location', 'northeast')
    end
    if i == 3; xlabel('Time (s)'); end

    % UKF (right column)
    subplot(3, 2, 2*i)
    hold on
    plot(tmeas, ukf_error(:, i), 'b', 'LineWidth', 0.8)
    plot(tmeas,  sigma_bounds_ukf(:, i), 'r--', 'LineWidth', 1)
    plot(tmeas, -sigma_bounds_ukf(:, i), 'r--', 'LineWidth', 1)
    hold off
    grid on
    if i == 1
        title('UKF')
        legend('Error', '3\sigma', 'Location', 'northeast')
    end
    if i == 3; xlabel('Time (s)'); end
end

sgtitle('Position Estimation Error with 3\sigma Bounds')
exportgraphics(gcf, 'Images/fig2_ekf_vs_ukf.png', 'Resolution', 300)


%% Figure 3: RSS Position Error — All Estimators Including Warm Start
figure('Position', [100 100 700 400])

% Batch/initial-state estimators
rss_bls   = sqrt(sum((measured_states(:,1:3) - bls_traj(:,1:3)).^2, 2));
rss_rls   = sqrt(sum((measured_states(:,1:3) - rls_traj(:,1:3)).^2, 2));
rss_glsdc = sqrt(sum((measured_states(:,1:3) - glsdc_traj(:,1:3)).^2, 2));

% Sequential estimators (cold start)
rss_ekf = sqrt(sum(ekf_error(:,1:3).^2, 2));
rss_ukf = sqrt(sum(ukf_error(:,1:3).^2, 2));

% Warm-started sequential estimators
rss_ekf_warm = sqrt(sum((measured_states_post(:,1:3) - xhat_ekf_warm(:,1:3)).^2, 2));
rss_ukf_warm = sqrt(sum((measured_states_post(:,1:3) - xhat_ukf_warm(:,1:3)).^2, 2));

semilogy(tmeas, rss_bls,   'c',   'LineWidth', 1.2, 'DisplayName', 'BLS')
hold on
semilogy(tmeas, rss_rls,   'm',   'LineWidth', 1.2, 'DisplayName', 'RLS-FF')
semilogy(tmeas, rss_glsdc, 'g',   'LineWidth', 1.2, 'DisplayName', 'GLSDC')
semilogy(tmeas, rss_ekf,   'b',   'LineWidth', 1.2, 'DisplayName', 'EKF')
semilogy(tmeas, rss_ukf,   'r',   'LineWidth', 1.2, 'DisplayName', 'UKF')
semilogy(tmeas_post, rss_ekf_warm, 'b--', 'LineWidth', 1.2, 'DisplayName', 'EKF (warm)')
semilogy(tmeas_post, rss_ukf_warm, 'r--', 'LineWidth', 1.2, 'DisplayName', 'UKF (warm)')
hold off

xlabel('Time (s)')
ylabel('RSS Position Error (km)')
title('Position Estimation Convergence — All Methods')
legend('Location', 'northeast')
grid on

exportgraphics(gcf, 'Images/fig3_rss_convergence.png', 'Resolution', 300)


%% Figure 4: RMSE Bar Chart — All Methods + Cold vs Warm
figure('Position', [100 100 800 350])

% --- Left: all methods, full arc ---
subplot(1, 2, 1)
names_all = categorical({'BLS', 'RLS-FF', 'GLSDC', 'EKF', 'UKF'});
names_all = reordercats(names_all, {'BLS', 'RLS-FF', 'GLSDC', 'EKF', 'UKF'});
pos_all   = [bls_pos; rls_pos; glsdc_pos; ekf_pos; ukf_pos];
vel_all   = [bls_vel; rls_vel; glsdc_vel; ekf_vel; ukf_vel];

bar(names_all, [pos_all, vel_all])
ylabel('RMSE')
title('Full Arc')
legend('Position (km)', 'Velocity (km/s)', 'Location', 'northwest')
grid on

% --- Right: cold vs warm, tail window (apples to apples) ---
subplot(1, 2, 2)
names_warm = categorical({'EKF (cold)', 'EKF (warm)', 'UKF (cold)', 'UKF (warm)'});
names_warm = reordercats(names_warm, {'EKF (cold)', 'EKF (warm)', 'UKF (cold)', 'UKF (warm)'});
pos_warm   = [ekf_tail_pos; ekf_warm_pos; ukf_tail_pos; ukf_warm_pos];
vel_warm   = [ekf_tail_vel; ekf_warm_vel; ukf_tail_vel; ukf_warm_vel];

bar(names_warm, [pos_warm, vel_warm])
ylabel('RMSE')
title('Tail Window — Cold vs Warm')
legend('Position (km)', 'Velocity (km/s)', 'Location', 'northwest')
grid on

sgtitle('Trajectory RMSE Comparison')
exportgraphics(gcf, 'Images/fig4_rmse_comparison.png', 'Resolution', 300)