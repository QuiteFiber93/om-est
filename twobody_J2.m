function ydot = twobody_J2(t, y, mu)
% Two-body dynamics with J2 oblateness perturbation, for use as a TRUTH model.
% State y = [r; v] (6x1), units km and km/s. ECI frame.
%
% J2 is the dominant gravitational perturbation for low Earth orbit. Using
% this as truth while the filters assume pure two-body creates a deliberate
% dynamics mismatch -- the process noise Q must absorb the unmodeled J2
% acceleration (see notes below).

    % Earth constants (consistent with mu = 398600.4415 km^3/s^2)
    J2   = 1.08262668e-3;     % Earth's second zonal harmonic (dimensionless)
    R_E  = 6378.137;          % Earth equatorial radius, km

    r = y(1:3);
    v = y(4:6);
    r_mag = norm(r);

    x = r(1); yy = r(2); z = r(3);

    % Point-mass (two-body) acceleration
    a_twobody = -mu / r_mag^3 * r;

    % J2 perturbing acceleration (ECI components)
    % Standard form: a_J2 = -(3/2) J2 (mu/r^2) (R_E/r)^2 * [ ... ]
    factor = -1.5 * J2 * (mu / r_mag^2) * (R_E / r_mag)^2;
    zr2 = (z / r_mag)^2;

    a_J2 = factor * [ (1 - 5*zr2) * x / r_mag;
                      (1 - 5*zr2) * yy / r_mag;
                      (3 - 5*zr2) * z / r_mag ];

    ydot = [v; a_twobody + a_J2];
end