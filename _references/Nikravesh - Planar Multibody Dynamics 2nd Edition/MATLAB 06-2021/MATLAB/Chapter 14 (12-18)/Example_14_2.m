% Filename: Example_14_2
clc; clear all

global s_O1_local u2_local uy theta0 k_r L0 k dc m J g

% Define body-fixed vectors, spring-damper, and inertia data
    s_O1_local = [0; 0.5]; u2_local = [0; -1];
    theta0 = 0; k_r = 6; 
    L0 = 0.7; k = 20; dc = 5; 
    m = [2; 1]; J = [0.04; 0.1]; g = 9.81; uy = [0; 1];
% Define joint coordinates and velocities
    theta = [pi/3; 0.8];
    theta_d = [0; 0];
% Integration array
    u = [theta; theta_d];
% Time span
    Tspan = 0:0.02:5.0;
% 4th order R-K
    [T, uT] = RK4(@Ex_14_2, Tspan, u);
    
    
subplot(1,2,1)
plot(T,uT(:,1))
hold on
subplot(1,2,2)
plot(T,uT(:,2))
hold on

    