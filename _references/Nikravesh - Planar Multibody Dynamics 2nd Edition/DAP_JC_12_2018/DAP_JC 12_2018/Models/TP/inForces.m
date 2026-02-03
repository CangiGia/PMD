function inForces
    include_global

F2 = Force_struct;
F2.type = 'weight';  % include the weight
F2.gravity = 9.81;   % default
F2.wgt = [0; -1]; % default

Forces = [F2];