% Filename: Example_13_3
clc; clear all

global s_O1_local u2_local uy theta0 k_r L0 k dc m M g

% Define body-fixed vectors, spring-damper, and inertia data
    s_O1_local = [0; 0.5]; u2_local = [0; -1]; uy = [0; 1];
    theta0 = 0; k_r = 6; 
    L0 = 0.7; k = 20; dc = 5; 
    m = [2; 1]; J = [0.04; 0.1]; g = 9.81;
% Mass matrix (array) (body coordinates)
    M = [m(1) m(1) J(1) m(2) m(2) J(2)]';
% Initial conditions for joint coordinates and velocities
    theta = [pi/3; 0.8];
    theta_d = [0; 0];
% Integration array
    u = [theta; theta_d];
% Time span
    Tspan = 0:0.02:5.0;
% 4th order R-K
    [T, uT] = RK4(@Ex_13_3, Tspan, u);
    