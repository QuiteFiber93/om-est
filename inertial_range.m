function rho = inertial_range(r, R, obsv_lat, LST)
    transformation = R * [cos(obsv_lat)*cos(LST); cos(obsv_lat)*sin(LST); sin(obsv_lat)*ones(size(LST))]';
    rho = r - transformation;
end