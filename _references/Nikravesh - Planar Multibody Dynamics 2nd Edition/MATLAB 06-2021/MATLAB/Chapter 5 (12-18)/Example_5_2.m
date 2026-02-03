% Example 5.2 (four-bar mechaniem)
    clc; clear all

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
    C = [L2*s_rot(u2) -L3*s_rot(u3)];
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
    Lpa = 0.2; beta2 = pi/6;
    upa = [cos(theta2 + beta2); sin(theta2 + beta2)];
    rP = L1*u1 + Lpa*upa 
    rP_d = L1*theta1_d*s_rot(u1) + Lpa*theta2_d*s_rot(upa)
    rP_dd = L1*(theta1_dd*s_rot(u1) - theta1_d^2*u1) + ...
            Lpa*(theta2_dd*s_rot(upa) - theta2_d^2*upa)  
