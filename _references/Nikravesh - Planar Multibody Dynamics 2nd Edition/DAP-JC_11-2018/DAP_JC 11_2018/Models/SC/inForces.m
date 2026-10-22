function inForces
    include_global

F1 = Force_struct;
F1.type = 'ptp'; % default
F1.iPindex = 1;
F1.jPindex = 6;
F1.k = 20;
F1.L0 = 0.6;

F2 = Force_struct;
F2.type = 'weight';  % include the weight
F2.gravity = 9.81;   % default
F2.wgt = [0; -1]; % default

Forces = [F1; F2];