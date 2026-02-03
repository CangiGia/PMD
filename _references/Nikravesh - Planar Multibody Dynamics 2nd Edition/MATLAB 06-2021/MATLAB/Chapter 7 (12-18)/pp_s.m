function f_s = pp_s(d, k, L0)
% This function computes the force of a point-to-point spring 
%    along vector d
% Compute length of d and unit vector u
    L = norm(d);
    u = d/L;
% Compute the spring force
    f = k*(L - L0);
% Compute the vector force along u
    f_s =  f*u;
end
