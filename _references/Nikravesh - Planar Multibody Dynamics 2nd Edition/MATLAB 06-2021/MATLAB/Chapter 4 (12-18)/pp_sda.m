function f_sda = pp_sda(d, d_d, k, L0, dc, fa)
% Compute the force of a point-to-point combined 
%    spring-damper-actuator element along vector d
% Compute length of d and unit vector u
    L = norm(d);
    u = d/L;
% Compute L_d
    L_d = d'*d_d/L;
% Compute spring and damper forces
    fs = k*(L - L0);
    fd = dc*L_d;
% Compute the sum of forces
    f = fs + fd + fa;
% Compute the vector force along u
    f_sda =  f*u;
end
