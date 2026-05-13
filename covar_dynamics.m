function ydot = covar_dynamics(t, y, Q, mu)
r = y(1:3);
v = y(4:6);
xdot = [v; -mu/norm(r)^3*r];

% Computing State Jacobian
r_mag = norm(r);
F = zeros(6, 6);
% dv/dv = I
F(1:3, 4:6) = eye(3);
F(4:6, 1:3) = 3*mu/r_mag^5 * (r*r') - mu/r_mag^3 * eye(3);

% Covariance Dynamics
P = reshape(y(7:42), 6, 6);
Pdot = F * P + P * F' + Q;

% Combined state dynamics
ydot = [xdot; Pdot(:)];
end