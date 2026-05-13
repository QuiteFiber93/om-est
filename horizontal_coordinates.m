function h_cords = horizontal_coordinates(rho_obsv)
    rho_mag = vecnorm(rho_obsv, 2, 2);
    az = wrapTo2Pi(atan2(rho_obsv(:, 2), rho_obsv(:, 3)));
    el = asin(rho_obsv(:, 1)./rho_mag);
    h_cords = [rho_mag, az, el];
end