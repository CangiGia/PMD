function ud = Ex_14_3_a(t, u) 
% Example 14.3(a): Coveyor Belt and friction

    global  m k L0 fN mu_s mu_d v_t v_conv f_eval

% Recover coordinates and velocities
    c = u(1); c_d = u(2);
% Force of the point-to-point spring 
    f_s = k*(c - L0);
% Relative velocity
    v_ij = v_conv - c_d;
% B-M friction model w/o viscous damping
    f_f = Friction_B(mu_s, mu_d, 0, v_t, 0.1, v_ij, fN);
% Total force
    fx = -f_s + f_f;
% Equations of motion and solution
    c_dd = fx/m;
% Construct ud array
    ud = [c_d; c_dd];
% Increment number of function evaluations
    f_eval = f_eval + 1;
end
