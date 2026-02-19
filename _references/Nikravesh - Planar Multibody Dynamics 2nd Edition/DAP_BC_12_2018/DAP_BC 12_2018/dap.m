%
% DAP_BC - Dynamic Analysis Program with Body Coordinate Formulation
%             Last update: August 2018
%      Complementary to Chapters 7 and 8 of the textbook
%
%             PLANAR MULTIBODY DYNAMICS
% Formulation, Programming with MATLAB, and Application
%                  Second Edition


    clear all

    addpath(genpath('Formulation'))
    folder = input(' Which folder contains the model? ', 's');
    addpath(fullfile(pwd,'Models',folder))
    
    include_global
    inData          % Read data
    initialize      % Determine numbers; pointers; dimensions; etc.

if nConst ~= 0
    ans = input('Do you want to correct the initial conditions? [(y)es/(n)o] ', 's');
    if ans == 'y'
        ic_correct  % Correct initial conditions
    end
    D = Jacobian;
    redund = rank(D); % Determine redundancy between constraints
    if redund < nConst
        disp('Redundancy in the constraints  ')
    end
end

% Pack coordinates and velocities ito u array
%   (ri, pi); for i = 1:nB and (ri_dot, pi_dot); for i = 1:nB 
    u = zeros(nB6,1); Bodies_to_u;
    t_initial = 0; t_final = input('Final time = ? '); % Input time data
    showtime = 1;
    if t_final == 0
        analysis(0, u); % one step function evaluation
        T = 0; uT = u';
    else
        dt    = input('Reporting time-step = ? ');
        Tspan = [t_initial:dt:t_final]; 
        options = odeset('RelTol', 1e-6, 'AbsTol', 1e-9);
        [T, uT] = ode45('analysis', Tspan, u, options);        
    end
    
    disp(['Number of function evaluations = ' num2str(num)])
    beep

    

%% Saving results to allow comparison with Python 
% filename = [folder '.txt'];
% output_data = [T uT];

% writematrix(output_data, filename);