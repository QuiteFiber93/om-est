function rho_obsv = observer_range(rho, obsv_lat, LST)
    rho_obsv = zeros(size(rho));
    rotation2 = [cos(obsv_lat), 0, sin(obsv_lat); 0, 1, 0; -sin(obsv_lat), 0, cos(obsv_lat)];
    for k =1:length(LST)
        LST_k = LST(k);
        rotation1 = [cos(LST_k), sin(LST_k), 0; -sin(LST_k), cos(LST_k), 0; 0, 0, 1];
        rho_obsv(k, :) = (rotation2 * rotation1 * rho(k, :)')';
    end
end