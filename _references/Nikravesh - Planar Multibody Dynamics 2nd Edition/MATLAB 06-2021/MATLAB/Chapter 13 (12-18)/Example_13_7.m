% Filename: Example_13_7
clc; 
clear all

global s_O1_local u1_local u2_local uy theta0 k_r L0 k dc m J g
global mc zeta omega Phi Phi_d D gamma

% Define body-fixed vectors
    s_O1_local = [0; 0.5]; u1_local = [0; -1]; u2_local = [0; -1];
% Define spring-damper constants and inertia data
    theta0 = 0; k_r = 6; 
    L0 = 0.7; k = 20; dc = 5; 
    m = [2; 1]; J = [0.04; 0.1]; g = 9.81; uy = [0; 1];
% Parameters for the penalty method
    mc = 1000; zeta = 1; omega = 20*pi;
% Define initial conditions for coordinates and velocities
    c = [0.4330; -0.2500; pi/3; 0.6928; -0.4000; pi/3];
    c_d = [0; 0; 0; 0; 0; 0];
% Integration array
    u = [c; c_d];
% Time span
    Tspan = 0:0.02:5.0;
% 4th order R-K
    [T, uT] = ode45(@Ex_13_7, Tspan, u);
    
% Recover constraint violations and Lagrange multipliers
    nt = length(T);
    Phi_norm   = zeros(nt,1); Phi_d_norm = zeros(nt,1); Lam = zeros(nt,4);
    for i = 1:nt
        ud = Ex_13_7(T(i),uT(i,:)');
        Phi_norm(i) = norm(Phi); Phi_d_norm(i) = norm(Phi_d);
        Lam(i,:) = -mc*(D*ud(7:12) - gamma + 2*zeta*omega*Phi_d + omega^2*Phi);
    end
% Plot violations        
    subplot(1,2,1)
    plot(T,Phi_norm)
    hold on
    subplot(1,2,2)
    plot(T,Phi_d_norm)
    hold on
% Plot Lagrange multipliers        
    figure
    subplot(1,2,1)
    plot(T,Lam(:,1))
    hold on
    plot(T,Lam(:,2),'g')
    subplot(1,2,2)
    plot(T,Lam(:,3),'r')
    hold on
    plot(T,Lam(:,4),'k')
