function inForces
    include_global

S1 = Force_struct;
S1.type = 'ptp'; % front spring-damper
S1.iPindex = 15;
S1.jPindex = 3;
S1.k = 75000;
S1.L0 = 0.29;
S1.dc = 6700;
S1.f_a = 0;

S2 = Force_struct;
S2.type = 'ptp'; % rear spring-damper
S2.iPindex = 16;
S2.jPindex = 17;
S2.k = 120000;
S2.L0 = 0.25;
S2.dc = 8500;
S2.f_a = 0;

S3 = Force_struct;
S3.type = 'user'; % tires
S3.k = 175000;
S3.L0 = 0.31;
S3.dc = 11000;

S4 = Force_struct;
S4.type = 'weight'; % 
S4.gravity = 9.81; % default
S4.wgt = [0; -1]; % default

Forces = [S1; S2; S3; S4];
