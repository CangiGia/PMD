%
% DAP_JC - Dynamic Analysis Program with Joint Coordinate Formulation
%             Last update: August 2018
%      Complementary to Chapter 9 of the textbook
%
%             PLANAR MULRIBODY DYNAMICS
% Formulation, Programming with MATLAB, and Application
%                  Second Edition


    clear all

    addpath(genpath('Formulation'))
    folder = input(' Which folder contains the model? ', 's');
    addpath(fullfile(pwd,'Models',folder))
    
    include_global
    global theta theta_d
    inData          % Read data
    initialize      % Determine numbers; pointers; dimensions; etc.

if nConst ~= 0
    ans = input('Do you want to correct the initial conditions? [(y)es/(n)o] ', 's');
    if ans == 'y'
        ic_correct  % Correct initial conditions
    end
    CutJoint(1, theta, theta_d, 0);
    redund = rank(C); % Determine redundancy between constraints
    if redund < nConst
        disp('Redundancy in the constraints  ')
    end
end

% Pack joint coordinates and velocities ito u array
%   theta(i); for i = 1:nJC and theta_d(i); for i = 1:nJC 
    u = [theta; theta_d]; 
    t_initial = 0; % t_final = 0;
    t_final = input('Final time = ? '); % Input time data
    showtime = 1;
    tSim = tic;
    if t_final == 0
        u_d = analysis(0, u); % one step function evaluation
        T = 0; uT = u';
    else
        dt    = input('Reporting time-step = ? ');
        Tspan = [t_initial:dt:t_final]; 
        options = odeset('RelTol', 1e-6, 'AbsTol', 1e-9);
        [T, uT] = ode45('analysis', Tspan, u, options);        
    end
    sim_time = toc(tSim); % wall-clock duration of the simulation [s]

% Position-constraint violation, ||Phi||, at every reporting time
    phi_norm = zeros(length(T), 1);
    if nConst ~= 0
        for k = 1:length(T)
            th = uT(k,1:nJC)'; th_d = uT(k,nJC+1:nJC2)';
            ChainOpen(1, th, th_d, T(k)); CutJoint(1, th, th_d, T(k));
            phi_norm(k) = norm(Phi);
        end
    end
    phi_max = max(phi_norm);
    
    disp(['Number of function evaluations = ' num2str(num)])
    disp(['Simulation time = ' num2str(sim_time) ' s'])
    disp(['Max constraint violation = ' num2str(phi_max)])
    disp(['Final constraint violation = ' num2str(phi_norm(end))])
    beep
