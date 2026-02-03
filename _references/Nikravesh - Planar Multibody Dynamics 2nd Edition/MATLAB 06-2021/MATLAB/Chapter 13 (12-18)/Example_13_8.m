% Filename: Example_13_8
clc; clear all

global s_O1_local u1_local u2_local uy theta0 k_r L0 k dc m M g

% Define body-fixed vectors
    s_O1_local = [0; 0.5]; u1_local = [0; -1]; u2_local = [0; -1];
% Define spring-damper constants and inertia data
    theta0 = 0; k_r = 6; 
    L0 = 0.7; k = 20; dc = 5; 
    m = [2; 1]; J = [0.04; 0.1]; g = 9.81; uy = [0; 1];
% Mass matrix (array) (body-coordinates)
    M = [m(1) m(1) J(1) m(2) m(2) J(2)]';
% Define initial conditions for coordinates and velocities
    c = [0.4330; -0.2500; pi/3; 0.6928; -0.4000; pi/3];
    c_d = zeros(6,1); mom = M.*c_d;
% Integration array
    u = [c; mom];
% Time span
    Tspan = 0:0.02:5.0;
% 4th order R-K (RK4) or ode45
    [T, uT] = ode45(@Ex_13_8, Tspan, u);
    