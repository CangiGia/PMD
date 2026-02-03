% Example 12.1 (slider-crank mechanism)
%     Example 5.1 (revised)
clc; clear all

% Constant data
    L1 = 0.12; L2 = 0.26;
% Driver data
    theta1_0 = 65*pi/180; omega1 = 2*pi; alpha1 = 0;
% Time parameter data
    del_t = 0.001; t_end = 1.0;
    nsteps = t_end/del_t + 1;
% arrays to save results
    T = zeros(nsteps,1);
    thetas    = zeros(nsteps,3);
    thetas_d  = zeros(nsteps,3);
    thetas_dd = zeros(nsteps,3);
for i = 1:nsteps 
    t = (i - 1)*del_t; T(i) = t;
    theta1 = theta1_0 + omega1*t + 0.5*alpha1*t^2;
% Solve for the other two variables (coordinates)
    theta2 = asin(-L1*sin(theta1)/L2);
    theta3 = L1*cos(theta1) + L2*cos(theta2);
    thetas(i,:) = [theta1 theta2 theta3];
% Unit vectors
    u1 = [cos(theta1); sin(theta1)]; u1r = s_rot(u1);
    u2 = [cos(theta2); sin(theta2)]; u2r = s_rot(u2);
% Jacibian matrix
    C = [L2*s_rot(u2)  [-1; 0]];
% Right-hand-side of velocity constraints
    theta1_d = omega1 + alpha1*t;
    rhsv = -L1*u1r*theta1_d;
% Solve for unknown velocities
    solution_v = C\rhsv;
    theta2_d = solution_v(1); theta3_d = solution_v(2);
    thetas_d(i,:) = [theta1_d theta2_d theta3_d];
% Right-hand-side of acceleration constraints
    theta1_dd = alpha1;
    rhsa = -L1*u1r*theta1_dd + L1*u1*theta1_d^2 + L2*u2*theta2_d^2;
% Solve for unknown accelerations
    solution_a = C\rhsa;
    theta2_dd = solution_a(1);
    theta3_dd = solution_a(2);
    thetas_dd(i,:) = [theta1_dd theta2_dd theta3_dd];
end
