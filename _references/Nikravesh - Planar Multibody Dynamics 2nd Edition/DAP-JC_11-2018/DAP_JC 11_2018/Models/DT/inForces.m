function inForces
    include_global

S1 = Force_struct;
S1.type = 'weight'; % 
S1.gravity = 9.81; % default
S1.wgt = [0; -1]; % default

Forces = [S1];
