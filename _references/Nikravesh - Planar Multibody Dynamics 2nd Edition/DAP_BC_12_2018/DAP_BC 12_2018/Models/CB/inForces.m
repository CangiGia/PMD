function inForces
    include_global

F1 = Force_struct;
F1.type = 'ptp'; % default
F1.iPindex = 2;
F1.jPindex = 1;
F1.k = 10;
F1.L0 = 0.8;
F1.dc = 0;

F2 = Force_struct;
F2.type = 'weight';  % include the weight

F3 = Force_struct;
F3.type = 'user';  

Forces = [F1; F2; F3];
