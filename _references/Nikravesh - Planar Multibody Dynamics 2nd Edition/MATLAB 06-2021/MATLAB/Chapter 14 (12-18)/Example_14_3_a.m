% Filename: Example_14_3_a
% Example 14.3(a): Coveyor Belt and friction
    clc; clear all

global  m k L0 fN mu_s mu_d v_t v_conv f_eval

    m = 1; g = 9.81; k = 10; L0 = 0.5; v_conv = 0.1;
    mu_s = 0.2; mu_d = 0.15; v_t = 0.001;
    fN = g*m;
% Define initial conditions
    c = 0.5; c_d = 0.1;
% Integration array
    u = [c; c_d];
% Time span
    Tspan = 0:0.02:10.0;
% Integrator ode45
    options = odeset('RelTol', 1e-6, 'AbsTol', 1e-9);
    f_eval = 0;
    [T, uT] = ode45(@Ex_14_3_a, Tspan, u, options);
    function_eval = f_eval
    subplot(2,1,1)
    plot(T,uT(:,1),'k')
    subplot(2,1,2)
    plot(T,uT(:,2),'k')
    