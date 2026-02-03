% Example 7.2 (two constrained bodies; continuation of Example 7.1)
    clc; clear all
    
%% Example 7.1
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
    f_s = pp_s(d1, k, L0);
    f_sd = pp_sd(d2, d2_d, 0, 0, dc);
    f_A1 = -f_s; f_C2 =  f_s; f_B1 = -f_sd;
% Compute moment of each force
    ns1 = s_rot(s_A1)'*f_A1;
    ns2 = s_rot(s_C2)'*f_C2;
    nd1 = s_rot(s_B1)'*f_B1;
% Construct vectors of the weights
    w1 = -m(1)*g*uy; w2 = -m(2)*g*uy;
% Construct array of applied forces
    h = [(f_A1 + f_B1 + w1)
         (ns1 + nd1)
         (f_C2 + w2)
         ns2];
    m_array = [m(1); m(1); J(1); m(2); m(2); J(2)];
    M = diag(m_array);
% Compute accelerations
%     acc = M\h

%% Example 7.2
% Additional data
    s_D2_local = [ 0; -0.12]; L = 0.355;
% Compute s_D2 and r_D2
    s_D2 = A2*s_D2_local; r_D2 = r2 + s_D2;
% Compute velocity of D2
    r_D2_d = r_Point_d(r2_d, s_D2, phi2_d);
% Construct vector d between B and D, its time derivative,
% and a unit vector
    d3   = r_B1 - r_D2; d3_d = r_B1_d - r_D2_d;
    u3 = d3/L;
% Evaluate the revolute-revolute joint constraint
    Phi = (u3'*d3 - L)/2
% Construct the Jacobian matrix
    D = [u3' u3'*s_rot(s_B1) -u3' -u3'*s_rot(s_D2)]
% Test for constraint violation at the velocity level
    c_d = [r1_d; phi1_d; r2_d; phi2_d];
    Phi_d = D*c_d
% Compute the time derivative of vector d and the unit vector
    d3_d = r_B1_d - r_D2_d; u3_d = d3_d/L;
% Compute the time derivative of s_B1 and s_D2
    s_B1_d = s_rot(s_B1)*phi1_d;
    s_D2_d = s_rot(s_D2)*phi2_d;
% Construct the r-h-s of acceleration constraint
    gamma = -u3_d'*d3_d + s_rot(u3)'* ...
            (s_B1_d*phi1_d - s_D2_d*phi2_d)
% Compute accelerations and Lagrange multiplier
    solution = [M -D'
                D  0 ]\[h; gamma]
