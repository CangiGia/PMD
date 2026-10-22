function inForces
    include_global

S1 = Force_struct;
S1.type = 'ptp'; % rear spring-damper
S1.iPindex = 11;
S1.jPindex = 12;
S1.k = 120000;
S1.L0 = 0.23;
S1.dc = 5300;
S1.f_a = 0;

S2 = Force_struct;
S2.type = 'ptp'; % front spring-damper
S2.iPindex = 14;
S2.jPindex = 13;
S2.k = 10000;
S2.L0 = 0.35;
S2.dc = 1350;
S2.f_a = 0;

S3 = Force_struct;
S3.type = 'user'; % tires
S3.k = 60000;
S3.L0 = 0.29;
S3.dc = 2700;

S4 = Force_struct;
S4.type = 'weight'; % 
S4.gravity = 9.81; % default
S4.wgt = [0; -1]; % default

Forces = [S1; S2; S3; S4];
