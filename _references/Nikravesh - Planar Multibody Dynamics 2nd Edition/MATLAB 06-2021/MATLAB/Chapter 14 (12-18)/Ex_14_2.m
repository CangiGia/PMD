function ud = Ex_14_2(t, u) 
% Example 14.2 equations of motion
global s_O1_local u2_local uy theta0 k_r L0 k dc m J g

% Recover joint coordinates and velocities
    theta = u(1:2); theta_d = u(3:4);
% Recursive coordinate transformation formulas
    phi1 = theta(1); 
    A1   = A_matrix(phi1); s_O1 = A1*s_O1_local;
    r1   = -s_O1;
    phi2 = phi1;
    A2   = A1;
    u2   = A2*u2_local;
    r2   = theta(2)*u2;
% Construct the B matrix
    z21 = [0; 0];
    B = [s_rot(r1) z21
         1          0
         s_rot(r2) u2
         1          0];
% Compute body velocities
    c_d = B*theta_d;
% Compute B_dot*theta_dot
    r1_d = c_d(1:2); r2_d = c_d(4:5);
    u2_d = s_rot(u2)*theta_d(1);
    Bdtd = [s_rot(r1_d)*theta_d(1)
            0
            s_rot(r2_d)*theta_d(1) + u2_d*theta_d(2)
            0];
% Construct body-coordinate's mass matrix
    M = [m(1) m(1) J(1) m(2) m(2) J(2)]';
% Compute the torque of the rotational spring 
    T_r = r_s(theta(1), k_r, theta0);
    T1 = -T_r;
% Compute the force of the point-to-point spring-damper 
    f_sd = pp_sd(r2, r2_d, k, L0, dc);
    f2 = -f_sd;
% Compute the gravitational forces
    w1 = -m(1)*g*uy; w2 = -m(2)*g*uy;
% Construct the array of applied forces (body coordinate)
    f_a = [w1; T1; w2 + f2; 0];
% Compute joint coordinate mass matrix and array of applied forces
    M_joint = B'*diag(M)*B;
    h_joint = B'*(f_a - diag(M)*Bdtd) - 10*theta_d; % fictitious damping
% determine joint accelerations
    theta_dd = M_joint\h_joint;
% Construct ud array
    ud = [theta_d; theta_dd];
end
