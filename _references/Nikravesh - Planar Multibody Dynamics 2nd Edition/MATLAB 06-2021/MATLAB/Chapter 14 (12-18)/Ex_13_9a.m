function ud = Ex_13_9a(t, u) 

global  L0 k dc m J weight f_eval
global  mu_d mu_s k_v v_s p k_t

    f_eval = f_eval + 1;
% Recover  coordinates and velocities
    c = u(1); c_d = u(2);
% Anderson et al. friction model
    v = c_d;
if abs(v) > v_s
    ff = Friction_A(mu_s, mu_d, v_s, p, k_t, v);
    fx = 8*sin(pi*t) - weight*(ff + k_v*v);
    c_dd = fx/m;
else
    fx = 8*sin(pi*t);
    DMD = [m -1
           1  0];
    rhs = [fx; 0];
    solution = DMD\rhs;
    if abs(solution(2)) < weight*mu_s
        c_dd = solution(1);
    else
        c_dd = (fx - weight*mu_s*sign(v))/m;
    end
end

% Construct ud array
    ud = [c_d; c_dd];
    
    