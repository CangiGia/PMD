function f_sd = pp_sd(d, d_d, k, L0, dc)
% This function computes the force of a point-to-point combined 
%    spring-damper element along vector d
% Compute length of d and unit vector u
    L = norm(d);
    u = d/L;
% Compute L_d
    L_d = u'*d_d;
% Compute spring and damper forces
    f = k*(L - L0) + dc*L_d;
% Compute the vector force along u
    f_sd =  f*u;
end
