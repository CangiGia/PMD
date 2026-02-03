% Example 5.1 (slider-crank mechanism)
clc; clear all

% Constant data
    L1 = 0.12; L2 = 0.26;
    
% Crank angle
    theta1 = 65*pi/180;
% Solve for the other two variables (coordinates)
    theta2 = asin(-L1*sin(theta1)/L2)
%     theta2_deg = theta2*180/pi
    theta3 = L1*cos(theta1) + L2*cos(theta2)
% ----------------
% For the second solution use the following commented data
%     theta2 = pi- theta2
% %     theta2_deg = theta2*180/pi
%     theta3 = L1*cos(theta1) + L2*cos(theta2)
% ----------------

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
