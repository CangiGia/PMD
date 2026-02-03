% Example 9.5-b
clc; clear all

%% Example 9.5
% Define constants
    r_O = [0.15; 0.33]; uy = [0; 1]; 
    L1 = 0.208; L2 = 0.6; L3 = 0.195;
% Define body-fixed vectors
    s_O1_local = [-L1/2; 0]; s_B1_local = [L1/2; 0]; 
    s_Q2_local = [-L2/2; 0]; s_B3_local = [0; L3];
    u3_local = [1; 0]; 
% Define spring-damper constants and inertia data
    m = [0.2; 1.0; 0.5]; J = [0.004; 0.1; 0.05]; g = 9.81;
% Provide initial conditions
    theta = [345*pi/180; 0.2149; 0.4017];
    theta_d = [pi; 1.4453; 0.5816];
% Recursive coordinate transformation formulas
    phi1 = theta(1); A1 = A_matrix(phi1);
    s_O1 = A1*s_O1_local;
    r1 = r_O - s_O1;
    phi2 = theta(2);
    A2 = A_matrix(phi2);
    r2 = -A2*s_Q2_local;
    phi3 = phi2;
    A3 = A2;
    u3 = A3*u3_local;
    r3 = theta(3)*u3;
    c = [r1; phi1; r2; phi2; r3; phi3];
% Construct the B matrix
    z21 = [0; 0];
    B = [s_rot(-s_O1)   z21         z21
         1              0           0
         z21            s_rot(r2) 	z21
         0              1           0
         z21            s_rot(r3) 	u3
         0              1           0]
% Compute body velocities
    c_d = B*theta_d;
% Compute B_dot*theta_dot
    r1_d = c_d(1:2); r2_d = c_d(4:5); r3_d = c_d(7:8);
    u3_d = s_rot(u3)*c_d(9);
    Bdtd = [s_rot(r1_d)*theta_d(1)
            0
            s_rot(r2_d)*theta_d(2)
            0
            s_rot(r3_d)*theta_d(2) + u3_d*theta_d(3)
            0]
% Construct body-coordinate's mass matrix
    M = [m(1) m(1) J(1) m(2) m(2) J(2) m(3) m(3) J(3)]';
% Compute the gravitational forces
    w1 = -m(1)*g*uy; w2 = -m(2)*g*uy; w3 = -m(3)*g*uy;
% Construct the array of applied forces (body coordinate)
    h_a = [w1; 0; w2; 0; w3; 0];
% Compute joint coordinate mass matrix and array of applied forces
    M_joint = B'*diag(M)*B
    h_joint = B'*(h_a - diag(M)*Bdtd)
    
%% Example 9.5-b
% Construct the joint coordinate Jacobian
    u1 = [cos(theta(1)); sin(theta(1))]; u2 = [cos(theta(2)); sin(theta(2))];
    u3r = s_rot(u3);
    C = [L1*s_rot(u1) (L3*u3 - theta(3)*u3r) -u3]
% Determine C_dot*theta_dot
    Cdtd = -L1*u1*theta_d(1)^2 + (L3*u3r + theta(3)*u3)*theta_d(2)^2 - ...
            2*u3r*theta_d(2)*theta_d(3)
% determine joint accelerations and Lagrange multipliers
    CMC = [M_joint -C'
           C       zeros(2)]
    rhs = [h_joint; -Cdtd]
    sol = CMC\rhs;
    theta_dd = sol(1:3)
    LagMult = sol(4:5)
    