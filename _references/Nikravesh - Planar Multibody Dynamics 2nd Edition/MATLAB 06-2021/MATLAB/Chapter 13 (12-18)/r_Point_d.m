function r_P_d = r_Point_d(r_d, s_P, phi_d)
% This function computes the velocity of point P
    r_P_d = r_d + s_rot(s_P)*phi_d;
end
