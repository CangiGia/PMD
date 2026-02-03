function inForces
    include_global

S1 = Force_struct;
S1.type = 'ptp'; % default
S1.iPindex = 2; % B1
S1.jPindex = 4; % O0
S1.k = 20000;
S1.L0 = 0.34;
S1.dc = 1100;

S2 = Force_struct;
S2.type = 'user'; % tire
S2.k = 100000;
S2.L0 = 0.30;
S2.dc = 1000;

S3 = Force_struct;
S3.type = 'weight'; % include the weight

Forces = [S1; S2; S3];