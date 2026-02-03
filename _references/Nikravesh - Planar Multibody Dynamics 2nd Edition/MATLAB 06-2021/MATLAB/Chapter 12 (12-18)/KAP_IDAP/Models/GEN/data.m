% ===== user supplied data for a four-bar mechanism ======
% ---- General fourbar ----
% Link lengths
    RQ = [5;-1.0];  % L_0(x) and L_0(y) 
%     L = [2 6 4];   % Grashof
    L = [2 7 3];    % Non-Grashof
% Secondary point (if none, set LP(i) = 0 and beta(i) = 0)
    LP = [1.0  5.5  6.5];
    beta = [20  22.5  -30];          % in degrees
% Initial angles 
    theta = [120 30 60];   % in degrees
% Driver data: theta(1) is already defined
    omega1 = 1.0;          % in rad/sec
    
% Parameters for the animation window
% axis_limits = [x_min x_max y_min y_max]
    axis_limits = [-3 12 -5 10]; 

% Mass and moment of inertia
    m = [1 1 1]; J = [1 1 1];
% Position of mass centers
    LG = [1.5 2.5 2.0]; gamma = [10 10 15];
% Gravity
    g = 9.81; ug = [0; -1];
    