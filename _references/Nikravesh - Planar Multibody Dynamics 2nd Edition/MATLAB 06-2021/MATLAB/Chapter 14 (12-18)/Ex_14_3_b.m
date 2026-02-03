function ud = Ex_14_3_b(t, u) 
% Example 14.3(b): Coveyor Belt and revised friction

    global  m k L0 fN mu_s mu_d v_t v_conv f_eval

% Recover  coordinates and velocities
    c = u(1); c_d = u(2);
% Force of the spring 
    f_s = k*(c - L0);
% Relative velocity
    v_ij = v_conv - c_d;
if abs(v_ij) >= v_t
    % B-M friction model w/o viscous damping
    f_f = Friction_B(mu_s, mu_d, 0, v_t, 0.1, v_ij, fN);
    fx = -f_s + f_f;
    c_dd = fx/m;
else
    fx = -f_s;
    DMD = [m -1
           1  0];
    rhs = [fx; 0];
    solution = DMD\rhs;
    if abs(solution(2)) < fN*mu_s
        c_dd = solution(1);
    else
        c_dd = (fx + fN*mu_s*sign(solution(2)))/m;
    end
end
% Construct ud array
    ud = [c_d; c_dd];
% Increment number of function evaluations
    f_eval = f_eval + 1;
