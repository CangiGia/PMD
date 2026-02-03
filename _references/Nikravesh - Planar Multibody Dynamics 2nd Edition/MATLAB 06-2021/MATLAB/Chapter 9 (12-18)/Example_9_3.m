% Example 9.3
clc; clear all

% Define body-fixed vectors
    s_O1_local = [0; 0.5]; u2_local = [0; -1]; uy = [0; 1];
% Define spring-damper constants and inertia data
    theta0 = 0; k_r = 6; 
    L0 = 0.7; k = 20; dc = 5; 
    m = [2; 1]; J = [0.04; 0.1]; g = 9.81;
% Define joint coordinates (inertial) and velocities
    theta = [pi/3; 0.8];
    theta_d = [0; 0];
% Recursive coordinate transformation formulas
    phi1 = theta(1); 
    A1 = A_matrix(phi1);
    r1 = -A1*s_O1_local;
    phi2 = phi1;
    A2 = A1;
    u2 = A2*u2_local;
    r2 = theta(2)*u2;
    c = [r1; phi1; r2; phi2]
% Construct the B matrix
    z21 = [0; 0];
    B = [s_rot(r1) z21
         1          0
         s_rot(r2) u2
         1          0]
% Compute body velocities
    c_d = B*theta_d
% Compute B_dot*theta_dot
    r1_d = c_d(1:2); r2_d = c_d(4:5);
    u2_d = s_rot(u2)*c_d(6);
    Bdtd = [s_rot(r1_d)*theta_d(1)
            0
            s_rot(r2_d)*theta_d(1) + s_rot(u2_d)*theta_d(2)
            0]
% Construct body-coordinate's mass matrix
    M = [m(1) m(1) J(1) m(2) m(2) J(2)]';
% Compute the torque of the rotational spring 
    T_r = r_s(theta(1), k_r, theta0)
    T1 = -T_r;
% Compute the force of the point-to-point spring-damper 
    f_sd = pp_sd(r2, r2_d, k, L0, dc)
    f2 = -f_sd;
% Compute the gravitational forces
    w1 = -m(1)*g*uy; w2 = -m(2)*g*uy;
% Construct the array of applied forces (body coordinate)
    h_a = [w1; T1; w2 + f2; 0]
% Compute joint coordinate mass matrix and array of applied forces
    M_joint = B'*diag(M)*B
    h_joint = B'*(h_a - diag(M)*Bdtd)
% determine joint accelerations
    theta_dd = M_joint\h_joint
    