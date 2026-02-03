% Example 6.6 (slider-crank mechanism; continuation of Examples 5.1 and 6.4)
    clc; clear all

%% Example 5.1
% Constant data
    L1 = 0.12; L2 = 0.26;
% Crank angle
    theta1 = 65*pi/180;
% Solve for the other two variables (coordinates)
    theta2 = asin(-L1*sin(theta1)/L2)
%     theta2_deg = theta2*180/pi
    theta3 = L1*cos(theta1) + L2*cos(theta2)
%     theta2 = pi- theta2
% %     theta2_deg = theta2*180/pi
%     theta3 = L1*cos(theta1) + L2*cos(theta2)
% Crank angular velocity
    theta1_d = 1.6;
% Unit vectors
    u1 = [cos(theta1); sin(theta1)]; u1r = s_rot(u1);
    u2 = [cos(theta2); sin(theta2)]; u2r = s_rot(u2);
% Jacibian matrix
    C = [L2*s_rot(u2)  [-1; 0]];
% Right-hand-side of velocity constraints
    rhsv = -L1*u1r*theta1_d;
% Solve for unknown velocities
    solution_v = C\rhsv;
    theta2_d = solution_v(1)
    theta3_d = solution_v(2)
    
% Crank angular acceleration
    theta1_dd = 0;
% Right-hand-side of acceleration constraints
    rhsa = -L1*u1r*theta1_dd + L1*u1*theta1_d^2 ...
                             + L2*u2*theta2_d^2;
% Solve for unknown accelerations
    solution_a = C\rhsa;
    theta2_dd = solution_a(1)
    theta3_dd = solution_a(2)

%% Example 6.4
% Constant data
    m1 = 1.0; J1 = 0.001; m2 = 2.0; J2 = 0.01; m3 = 4.0; J3 = 0.02;
    g = 9.81; uy = [0; 1]; fB = [50; 0];
% Gravitational forces
    w1 = -g*m1*uy; w2 = -g*m2*uy; w3 = -g*m3*uy;
% Array of "known" applied forces
    h_a = [w1; 0;w2; 0; (w3 + fB); 0]
% Position vectors (moment arms)
    sO1 = -(L1/2)*u1; sO1r = s_rot(sO1);
    sA1 =  (L1/2)*u1; sA1r = s_rot(sA1);
    sA2 = -(L2/2)*u2; sA2r = s_rot(sA2);
    sB2 =  (L2/2)*u2; sB2r = s_rot(sB2);
% Coefficient matrix for reaction forces
    I2 = eye(2); Z2 = zeros(2); Z12 = zeros(1,2); Z21 = zeros(2,1);
    Dt = [I2    -I2     Z2    Z21  Z21
          sO1r' -sA1r'  Z12   0    0
          Z2     I2    -I2    Z21  Z21
          Z12    sA2r' -sB2r' 0    0
          Z2     Z2     I2    uy   Z21
          Z12    Z12    Z12   0    1]
% Mass matrix
    M_diag = [m1 m1 J1 m2 m2 J2 m3 m3 J3]';
  
%% Example 6.6
% Coordinates of points
    rG1 = (L1/2)*u1; rA  = L1*u1; rG2 = rA + (L2/2)*u2; rB  = [theta3; 0];
% Acceleration of mass centers
    rG1_dd = s_rot(rG1)*theta1_dd - rG1*theta1_d^2;
    rG2_dd = s_rot(rA)*theta1_dd - rA*theta1_d^2 - sA2r*theta2_dd + sA2*theta2_d^2;
    rG3_dd = [theta3_dd; 0];
    c_dd = [rG1_dd; theta1_dd; rG2_dd; theta2_dd; rG3_dd; theta3_dd];
% Mass times acceleration
    Mcdd = M_diag.*c_dd
% Add a 9th column to Dt for the unknown torque
    column9 = zeros(9,1); column9(3) = 1.0;
    Dt9 = [Dt column9];
% Solve the equilibrium equations
    sol = Dt9\(Mcdd - h_a);
    T_a = sol(9)
