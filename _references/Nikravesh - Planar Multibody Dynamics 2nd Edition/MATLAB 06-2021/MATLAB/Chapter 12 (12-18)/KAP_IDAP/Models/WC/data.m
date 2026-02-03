% ===== user supplied data for a four-bar mechanism ======
% ---- Web_cutter ----

% Link lengths
    RQ = [0.65; -0.1];  % L_0(x) and L_0(y) 
    L = [0.1 0.7 1.0];     
% Secondary point (if none, set LP(i) = 0 and beta(i) = 0)
    LP = [0 1.0 0.82];
    beta = [0 -36.87 -37.57];          % in degrees
% Initial angles
    theta = [30 90 140];
% Driver data: theta(1) is already defined
    omega1 = 2*pi;          % in rad/sec
    
% Parameters for the animation window
% axis_limits = [x_min x_max y_min y_max]
    axis_limits = [-0.2 1 -0.2 1]; 

% Mass and moment of inertia
    m = [1 10 10]; J = [1.0 5.0 5.0];
% Position of mass centers
    LG = [0 0.52 0.68]; gamma = [0 -16.7 -17.1];
% Gravity
    g = 9.81; ug = [0; -1];