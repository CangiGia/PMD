function inForces
    include_global

S1 = Force_struct;
S1.type = 'rot-sda'; % neck pin 0
S1.iBindex = 2;
S1.jBindex = 1;
S1.k = 271;
S1.theta0 = 0;
S1.dc = 0.11;
S1.T_a = 0;

S2 = Force_struct;
S2.type = 'rot-sda'; % neck pin 1
S2.iBindex = 3;
S2.jBindex = 2;
S2.k = 271;
S2.theta0 = 0;
S2.dc = 0.11;
S2.T_a = 0;

S3 = Force_struct;
S3.type = 'rot-sda'; % neck pin 2
S3.iBindex = 4;
S3.jBindex = 3;
S3.k = 271;
S3.theta0 = 0;
S3.dc = 0.11;
S3.T_a = 0;

S4 = Force_struct;
S4.type = 'rot-sda'; % neck pin 3
S4.iBindex = 5;
S4.jBindex = 4;
S4.k = 271;
S4.theta0 = 0;
S4.dc = 0.11;
S4.T_a = 0;

S5 = Force_struct;
S5.type = 'rot-sda'; % neck pin 4
S5.iBindex = 6;
S5.jBindex = 5;
S5.k = 271;
S5.theta0 = 0;
S5.dc = 0.11;
S5.T_a = 0;

S6 = Force_struct;
S6.type = 'rot-sda'; % neck pin 5
S6.iBindex = 7;
S6.jBindex = 6;
S6.k = 271;
S6.theta0 = 0;
S6.dc = 0.11;
S6.T_a = 0;

S7 = Force_struct;
S7.type = 'rot-sda'; % neck pin 6
S7.iBindex = 8;
S7.jBindex = 7;
S7.k = 271;
S7.theta0 = 0;
S7.dc = 0.11;
S7.T_a = 0;

S8 = Force_struct;
S8.type = 'rot-sda'; % neck pin 7
S8.iBindex = 9;
S8.jBindex = 8;
S8.k = 271;
S8.theta0 = 0;
S8.dc = 0.22;
S8.T_a = 0;

S9 = Force_struct;
S9.type = 'rot-sda'; % brain-skull rotation
S9.iBindex = 10;
S9.jBindex = 9;
S9.k = 23;
S9.theta0 = 0;
S9.dc = 0.05;
S9.T_a = 0;

S10 = Force_struct;
S10.type = 'ptp'; % brain-skull
S10.iPindex = 19;
S10.jPindex = 23;
S10.k = 8750;
S10.L0 = 0.04;
S10.dc = 90;
S10.f_a = 0;

S11 = Force_struct;
S11.type = 'ptp'; % brain-skull
S11.iPindex = 20;
S11.jPindex = 23;
S11.k = 8750;
S11.L0 = 0.04;
S11.dc = 90;
S11.f_a = 0;

S12 = Force_struct;
S12.type = 'ptp'; % brain-skull
S12.iPindex = 21;
S12.jPindex = 23;
S12.k = 8750;
S12.L0 = 0.04;
S12.dc = 90;
S12.f_a = 0;

S13 = Force_struct;
S13.type = 'ptp'; % brain-skull
S13.iPindex = 22;
S13.jPindex = 23;
S13.k = 8750;
S13.L0 = 0.04;
S13.dc = 90;
S13.f_a = 0;

S14 = Force_struct;
S14.type = 'user'; % torso shaking force

S15 = Force_struct;
S15.type = 'weight'; % 
S15.gravity = 9.81; % default
S15.wgt = [0; -1]; % default

Forces = [S1; S2; S3; S4; S5; S6; S7; S8; S9; S10; S11; S12; S13; S14; S15];
