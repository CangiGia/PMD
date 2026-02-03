% ===== user supplied data for a four-bar mechanism ======
% ---- Film-strip Advancer ----
% Link lengths
    RQ = [-3; 4.5];  % L_0(x) and L_0(y) 
    L = [1.5 5  3];     
% Secondary point (if none, set LP(i) = 0 and beta(i) = 0)
    LP = [0 10.5 0];
    beta = [0 2.73 0];          % in degrees
% Initial angles
    theta = [30 90 30];          % in degrees
% Driver data: theta(1) is already defined
    omega1 = -1.0;          % in rad/sec
    
% Parameters for the animation window
% axis_limits = [x_min x_max y_min y_max]
    axis_limits = [-4 3 -2 13];
    
% Mass and moment of inertia
    m = [1 1 1]; J = [1 1 1];
% Position of mass centers
    LG = [0 0 0]; gamma = [0 0 0];
% Gravity
    g = 9.81; ug = [0; -1];
