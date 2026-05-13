function H_k = compute_H(x, obsv_lat, LST, R_obsv)

% H_k = [H_11, 0];
% H_11 = [drho/dx; daz/dx; del/dx]

% Defining measurement values
r_inertial = x(1:3);
rho_inertial = inertial_range(r_inertial', R_obsv, obsv_lat, LST);
rho_obsv = observer_range(rho_inertial, obsv_lat, LST);

% Extract components
rho_u = rho_obsv(1);
rho_e = rho_obsv(2);
rho_n = rho_obsv(3);
rho_mag = norm(rho_obsv);

% Preallocating H11
H11 = zeros(3, 3);

% First Row
H11(1, 1) = (rho_u * cos(obsv_lat) * cos(LST) - rho_e * sin(LST) - rho_n * sin(obsv_lat) * cos(LST)) / rho_mag;
H11(1, 2) = (rho_u * cos(obsv_lat) * sin(LST) + rho_e * cos(LST) - rho_n * sin(obsv_lat) * sin(LST)) / rho_mag;
H11(1, 3) = (rho_u * sin(obsv_lat) + rho_n * cos(obsv_lat)) / rho_mag;

% Second Row
denom_az = rho_n^2 + rho_e^2;
H11(2, 1) = (rho_e * sin(obsv_lat) * cos(LST) - rho_n * sin(LST)) / denom_az;
H11(2, 2) = (rho_e * sin(obsv_lat) * sin(LST) + rho_n * cos(LST)) / denom_az;
H11(2, 3) = -rho_e * cos(obsv_lat) / denom_az;

% Third Row
denom_el = rho_mag * sqrt(rho_mag^2 - rho_u^2);
H11(3, 1) = (rho_mag * cos(obsv_lat) * cos(LST) - rho_u * H11(1, 1)) / denom_el;
H11(3, 2) = (rho_mag * cos(obsv_lat) * sin(LST) - rho_u * H11(1, 2)) / denom_el;
H11(3, 3) = (rho_mag * sin(obsv_lat) - rho_u * H11(1, 3)) / denom_el;

H_k = [H11, zeros(3, 3)];
end