function ydot = twobody_STM(t, y, mu)
    r = y(1:3);
    v = y(4:6);
    Phi = reshape(y(7:42), 6, 6);

    % Phidot = df/dx * Phi
    % A = df/dx = [dv/dr, dv/dv; dvdot/dr, dvdot/dv]
    %           = [0, I; A_21, 0]
    A = [zeros(3), eye(3); 3*mu/norm(r)^5*(r*r') - mu/norm(r)^3*eye(3), zeros(3)];
    Phidot = A * Phi;
    ydot = [v; -mu*r/norm(r)^3; Phidot(:)];
end