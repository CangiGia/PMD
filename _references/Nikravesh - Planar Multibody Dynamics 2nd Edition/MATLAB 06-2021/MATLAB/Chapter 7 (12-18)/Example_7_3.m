% Example 7.3 (Variable-length pendulum)
    clc; clear all

% Define body-fixed vectors, spring-damper, and inertia data (constants)
    s_O1_local = [0; 0.5]; u1_local = [0; 1]; u2_local = [0; 1];
    theta0 = 0; k_r = 6; 
    L0 = 0.7; k = 20; dc = 5; 
    m = [2; 1]; J = [0.04; 0.1]; g = 9.81; uy = [0; 1];
% Define coordinates and velocities (variables)
    r1 = [0.4330; -0.2500]; phi1 = pi/3;
    r2 = [0.6928; -0.4000]; phi2 = pi/3;
    r1_d = [ 0; 0]; phi1_d = 0;
    r2_d = [ 0; 0]; phi2_d = 0;

% Compute A matrices
    A1 = A_matrix(phi1); A2 = A_matrix(phi2);
% Compute global components of vectors
    s_O1 = A1*s_O1_local;
    u1 = A1*u1_local; u2 = A2*u2_local;
% Construct Jacobian D
    d = r2 - r1; d_d = r2_d - r1_d; u1r = s_rot(u1);
    D = [eye(2)    s_rot(s_O1) zeros(2,3)
        -u1r' -u1'*d       u1r'   0
         0 0  -1           0 0    1]
% Construct r-h-s of acc constraints
    gamma = [s_O1*phi1_d^2
             (u1r'*d*phi1_d + 2*u1'*d_d)*phi1_d
                          0]
% Construct body-coordinate's mass matrix
    M = [m(1) m(1) J(1) m(2) m(2) J(2)]';
% Compute the torque of the rotational spring 
    T_r = r_s(phi1, k_r, theta0)
    T1 = -T_r;
% Compute the force of the point-to-point spring-damper 
    f_sd = pp_sd(r2, r2_d, k, L0, dc)
    f2 = -f_sd;
% Compute the gravitational forces
    w1 = -m(1)*g*uy; w2 = -m(2)*g*uy;
% Construct the array of applied forces
    f_a = [w1; T1; w2 + f2; 0]
% Construct mass/Jacobian matrix and array of forces/gamma
    DMD = [diag(M) -D'
           D  zeros(4)]
    rhs = [f_a; gamma]
% Determine accelerations and Lagrange multipliers
    solution = DMD\rhs;
    c_dd = solution(1:6)
    Lag  = solution(7:10)
