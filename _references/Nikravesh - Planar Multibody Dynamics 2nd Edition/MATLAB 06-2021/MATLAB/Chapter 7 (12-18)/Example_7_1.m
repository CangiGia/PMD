% Example 7.1 (two unconstrained bodies)
    clc, clear all
    
% Define mass and moment of inertia for the bodies
    m = [0.2 0.15]; J = [0.03 0.02]; g = 9.81; uy = [0; 1];
% Define spring-damper-actiator constants
    k = 50; L0 = 0.2; dc = 20;    
% Define body-fixed vectors
    s_A1_local = [ 0.15; 0]; s_B1_local = [-0.15; 0]; 
    s_C2_local = [ 0; 0.12]; r_Q0 = [ 0; 0]; 
% Define coordinates (inertial) and velocities
    r1 = [ 0.02; 0.2]; phi1 = pi/4;
    r2 = [ 0.22; 0.1]; phi2 = 15*pi/180;
    r1_d = [-0.05; 0.1] ; phi1_d = -0.3;
    r2_d = [-0.09; 0.13]; phi2_d = 0.07;

% Compute A matrices
    A1 = A_matrix(phi1); A2 = A_matrix(phi2);
% Compute global components of s_A1, s_B1 and s_C2
    s_A1 = A1*s_A1_local; s_B1 = A1*s_B1_local; s_C2 = A2*s_C2_local;
% Compute global coordinates of A1, B1 and C2
    r_A1 = r1 + s_A1; r_B1 = r1 + s_B1; r_C2 = r2 + s_C2;
% Compute velocity of A1, B1 and C2
    r_B1_d = r_Point_d(r1_d, s_B1, phi1_d);
% Compute vector d between A and C and its time derivative
    d1   = r_A1 - r_C2; 
% Compute vector d between B and Q and its time derivative
    d2   = r_B1 - r_Q0; d2_d = r_B1_d;
% Compute force of the spring and the force of the damper
    f_s = pp_s(d1, k, L0)
    f_sd = pp_sd(d2, d2_d, 0, 0, dc)
    f_A1 = -f_s; f_C2 =  f_s; f_B1 = -f_sd;
% Compute moment of each force
    ns1 = s_rot(s_A1)'*f_A1 
    ns2 = s_rot(s_C2)'*f_C2
    nd1 = s_rot(s_B1)'*f_B1
% Construct vectors of the weights
    w1 = -m(1)*g*uy; w2 = -m(2)*g*uy;
% Construct array of applied forces
    h = [(f_A1 + f_B1 + w1)
         (ns1 + nd1)
         (f_C2 + w2)
         ns2]
    m_array = [m(1); m(1); J(1); m(2); m(2); J(2)];
    M = diag(m_array)
% Compute accelerations
    acc = M\h
