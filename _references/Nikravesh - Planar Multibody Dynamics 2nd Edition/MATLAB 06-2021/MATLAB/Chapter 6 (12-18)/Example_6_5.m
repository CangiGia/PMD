% Example 6.5 (four-bar mechanism; continuation of Example 5.2)
    clc; clear all

%% Example 5.2
    global L0 L1 L2 L3
% Constant data
    L0 = 0.2; L1 = 0.1; L2 = 0.3; L3 = 0.22;
% Crank angle
    theta1 = pi/4;
% estimate for the other two angles
    theta23 = [0.5; 1.3]; % estimates
%     theta23 = [4.0; 4.0]; % estimates
% fsolve's parameters
    options = optimset('display', 'off');
    angles = fsolve(@fourbar, theta23, options, theta1);
    angles_deg = angles*180/pi;
    theta2 = angles(1)
    theta3 = angles(2)
% Crank angular velocity
    theta1_d = 1.5;
% Unit vectors
    u1 = [cos(theta1); sin(theta1)];
    u2 = [cos(theta2); sin(theta2)];
    u3 = [cos(theta3); sin(theta3)];
% Jacibian matrix
    C = [L2*s_rot(u2)  L3*s_rot(u3)];
% Right-hand-side of velocity constraints
    rhsv = -L1*s_rot(u1)*theta1_d;
% Solve for unknown velocities
    solution_v = C\rhsv;
    theta2_d = solution_v(1)
    theta3_d = solution_v(2)

% Solve for unknown velocities
    theta1_dd = -1.2;
% Right-hand-side of acceleration constraints
    rhsa = -L1*s_rot(u1)*theta1_dd + L1*u1*theta1_d^2 + ...
            L2*u2*theta2_d^2 - L3*u3*theta3_d^2;
% Solve for unknown accelerations
    solution_a = C\rhsa;
    theta2_dd = solution_a(1)
    theta3_dd = solution_a(2)
% Coupler Point 
    LPA = 0.2; beta2 = pi/6;
    uPA = [cos(theta2 + beta2); sin(theta2 + beta2)];
    rP = L1*u1 + LPA*uPA 
    rP_d = L1*theta1_d*s_rot(u1) + LPA*theta2_d*s_rot(uPA)
    rP_dd = L1*(theta1_dd*s_rot(u1) - theta1_d^2*u1) + ...
            LPA*(theta2_dd*s_rot(uPA) - theta2_d^2*uPA)  
       
%% Example 6.5
% Constant data
    LG2A = 0.16; gamma2 = 12*pi/180;
    m1 = 1.0; J1 = 0.001; m2 = 6.0; J2 = 0.1; m3 = 3.0; J3 = 0.015;
    g = 9.81; dc = 125;
    LRO = 0.2; uy = [0; 1];
% Gravitational forces
    w1 = -g*m1*uy; w2 = -g*m2*uy; w3 = -g*m3*uy;
% Damping force
    rR = [0; LRO];
    d = rP - rR; L = norm(d); u = d/L;
    L_d = u'*rP_d 
    fd = dc*L_d
    fP2 = -fd*u
% Position vectors
    rA  = L1*u1; rB  = rA + L2*u2; rQ = [L0; 0];
    uG2A = [cos(theta2 + gamma2); sin(theta2 + gamma2)]; sG2 = LG2A*uG2A;
    rG1 = (L1/2)*u1; rG2 = rA + sG2; rG3 = rQ + (L3/2)*u3; 
% Moment arms
    sO1 = -(L1/2)*u1; sO1r = s_rot(sO1);
    sA1 =  (L1/2)*u1; sA1r = s_rot(sA1);
    sA2 =  rA - rG2;  sA2r = s_rot(sA2);
    sB2 =  rB - rG2;  sB2r = s_rot(sB2);
    sP2 =  rP - rG2;  sP2r = s_rot(sP2);
    sB3 =  (L3/2)*u3; sB3r = s_rot(sB3);
    sQ3 = -(L3/2)*u3; sQ3r = s_rot(sQ3);
% Array of "known" applied forces
    h_a = [w1; 0; (w2 + fP2); sP2r'*fP2; w3; 0]
% Coefficient matrix for reaction forces
    I2 = eye(2); Z2 = zeros(2); Z12 = zeros(1,2); Z21 = zeros(2,1);
    Dt = [I2    -I2     Z2     Z2
          sO1r' -sA1r'  Z12    Z12    
          Z2     I2    -I2     Z2
          Z12    sA2r' -sB2r'  Z12
          Z2     Z2     I2     I2
          Z12    Z12    sB3r'  sQ3r']
% Mass matrix
    M_diag = [m1 m1 J1 m2 m2 J2 m3 m3 J3]';
    M = diag(M_diag);
    
    