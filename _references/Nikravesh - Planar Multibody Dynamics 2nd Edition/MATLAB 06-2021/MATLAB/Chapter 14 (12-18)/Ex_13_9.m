function ud = Ex_13_9(t, u) 
% Example 13.8: Coveyor Belt and friction

global  L0 k dc m J weight f_eval
global  mu_d mu_s k_v v_s p k_t

    f_eval = f_eval + 1;
% Recover  coordinates and velocities
    c = u(1); c_d = u(2);
% % Anderson et al. friction model
    v = c_d;
    ff = Friction_A(mu_s, mu_d, v_s, p, k_t, v);
    fx = 8*sin(pi*t) - weight*(ff + k_v*v);
% Equations of motion and solution
    c_dd = fx/m;
% Construct ud array
    ud = [c_d; c_dd];
    
    