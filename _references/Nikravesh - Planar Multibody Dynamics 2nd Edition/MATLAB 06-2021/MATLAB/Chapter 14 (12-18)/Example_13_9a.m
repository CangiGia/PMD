% Filename: Example_13_9a
% Moving block and revised Anderson friction model
clc; clear all

global  L0 k dc m J weight f_eval
global  mu_d mu_s k_v v_s p k_t

    f_eval = 0;
% Define spring-damper constants and inertia data
    L0 = 1.0; k = 10; dc = 1; 
    m = 1; J = 1; weight = 9.81;
    mu_d = 0.5; mu_s = 0.7; k_v = 0.0; v_s = 0.001; p = 2; k_t = 10000;
% Define initial conditions for coordinates and velocities
    c = 0; c_d = 0; % Block
% Integration array
    u = [c; c_d];
% Time span
    Tspan = 0:0.02:6;
% 4th order R-K (RK4) or ode45
    options = odeset('RelTol', 1e-6, 'AbsTol', 1e-9);
    [T, uT] = ode45(@Ex_13_9a, Tspan, u, options);
    function_eval = f_eval
    plot(T,uT(:,1),'g')
    hold on
    plot(T,uT(:,2),'r')
    
    