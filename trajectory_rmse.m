function [rmse_comp, rmse_pos, rmse_vel] = trajectory_rmse(xhat_traj, true_traj)
    err = xhat_traj - true_traj;                       % N x 6
    rmse_comp = sqrt(mean(err.^2, 1));                 % 1 x 6
    pos_err_mag = sqrt(sum(err(:,1:3).^2, 2));         % N x 1
    vel_err_mag = sqrt(sum(err(:,4:6).^2, 2));         % N x 1
    rmse_pos = sqrt(mean(pos_err_mag.^2));             % scalar
    rmse_vel = sqrt(mean(vel_err_mag.^2));             % scalar
end