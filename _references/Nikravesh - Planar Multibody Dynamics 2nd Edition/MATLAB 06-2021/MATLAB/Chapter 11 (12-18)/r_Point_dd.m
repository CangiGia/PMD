function r_P_dd = r_Point_dd(r_dd, s_P, phi_d, phi_dd)
% This function computes the acceleration of point P
    r_P_dd = r_dd + s_rot(s_P)*phi_dd - s_P*phi_d^2;
end
