function user_force
    include_global

% Anderson et al. friction model
    mu_d = 0.15; mu_s = 0.2; mu_v = 0.0; v_s = 0.001; p = 2; k_t = 10000;
    fy = 9.81; % normal force
    v_conv = 0.1;
    v = v_conv - Bodies(1).r_d(1);
    ff = Friction_A(mu_s, mu_d, v_s, p, k_t, v, fy);
    fx = ff + mu_v*v*fy;
    fs = [fx; 0];
    Bodies(1).f = Bodies(1).f + fs;
end
