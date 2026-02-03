function inForces
    include_global

S1 = Force_struct;
S1.type = 'ptp'; % default
S1.iPindex = 9;
S1.jPindex = 10;
S1.f_a = 0;
S1.k = 90000;
S1.L0 = 0.23;
S1.dc = 1100;

S2 = Force_struct;
S2.type = 'user'; % tire
S2.k = 50000;
S2.L0 = 0.35;
S2.dc = 1000;

S3 = Force_struct;
S3.type = 'weight';  % include the weight
S3.gravity = 9.81;   % default
S3.wgt = [0; -1]; % default

Forces = [S1; S2; S3];
