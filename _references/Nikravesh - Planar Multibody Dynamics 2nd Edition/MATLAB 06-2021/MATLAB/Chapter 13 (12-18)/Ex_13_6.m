function ud = Ex_13_6(t, u) 
% Example 13.6: Equations of motion

global s_O1_local u1_local u2_local uy theta0 k_r L0 k dc m M g
global alpha beta Phi Phi_d Lag

% Recover  coordinates and velocities
    c = u(1:6); c_d = u(7:12);
% Extract coordinates and velocities
    r1 = c(1:2); phi1 = c(3); A1 = A_matrix(phi1);
    r2 = c(4:5); phi2 = c(6); A2 = A1; 
    s_O1 = A1*s_O1_local; r_O1 = r1 + s_O1;
    u1 = A1*u1_local; u2 = A2*u2_local;
    r1_d = c_d(1:2); phi1_d = c_d(3); 
    r2_d = c_d(4:5); phi2_d = c_d(6); 
% Jacobian D
    d = r2 - r1; d_d = r2_d - r1_d; u1r = s_rot(u1);
    D = [eye(2)    s_rot(s_O1) zeros(2,3)
         0 0 -1            0 0   1
        -u1r' -u1'*d       u1r'   0];
% Right-hand side of acc constraints
    gamma = [s_O1*phi1_d^2
             0
             (u1r'*d*phi1_d + 2*u1'*d_d)*phi1_d];
% Torque of the rotational spring 
    T_r = r_s(phi1, k_r, theta0);
    T1 = -T_r;
% Force of the point-to-point spring-damper 
    f_sd = pp_sd(r2, r2_d, k, L0, dc);
    f2 = -f_sd;
% Gravitational forces
    w1 = -m(1)*g*uy; w2 = -m(2)*g*uy;
% Array of applied forces (body coordinates)
    f_a = [w1; T1; w2 + f2; 0];
% Evaluate constraints at coordinate and velocity levels
    Phi = [r_O1; u1r'*u2; u1r'*d];
    Phi_d = D*c_d;
% Include stabilization terms in gamma array
    gamma = gamma - 2*alpha*Phi_d - beta^2*Phi;
% Mass/Jacobian matrices and array of applied forces/gamma
    DMD = [diag(M) -D'
           D  zeros(4)];
    rhs = [f_a; gamma];
% Compute accelerations and Lagrange multipliers
    solution = DMD\rhs;
% Extract accelerations
    c_dd = solution(1:6);
    Lag  = solution(7:10);
% Construct ud array
    ud = [c_d; c_dd];
end
    