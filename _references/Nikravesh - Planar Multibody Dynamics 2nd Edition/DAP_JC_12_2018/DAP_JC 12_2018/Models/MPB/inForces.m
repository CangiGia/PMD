function inForces
    include_global
    Force = Force_struct;

S1 = Force;
S1.type = 'ptp'; % default
S1.iPindex = 4; % B1
S1.jPindex = 6; % O0
S1.k = 22600;
S1.L0 = 0.34;
S1.dc = 1270;
S1.f_a = 0;

S2 = Force;
S2.type = 'user'; % tire
S2.k = 200000;
S2.L0 = 0.30;
S2.dc = 2000;

S3 = Force;
S3.type = 'weight';  % include the weight

Forces = [S1; S2; S3];