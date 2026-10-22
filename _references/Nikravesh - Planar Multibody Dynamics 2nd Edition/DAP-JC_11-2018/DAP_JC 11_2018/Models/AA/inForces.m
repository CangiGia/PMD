function inForces
    include_global
    Force = Force_struct;

S1 = Force;
S1.type = 'ptp'; % default
S1.iPindex = 9;
S1.jPindex = 10;
S1.k = 91600;
S1.L0 = 0.23;
S1.dc = 1433;
S1.f_a = 0;

S2 = Force;
S2.type = 'user'; % tire
S2.k = 58770;
S2.L0 = 0.35;
S2.dc = 2500;

S3 = Force;
S3.type = 'weight';  % include the weight
S3.gravity = 9.81;   % default
S3.wgt = [0; -1]; % default

Forces = [S1; S2; S3];
