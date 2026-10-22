function inPoints
    include_global

P1 = Point_struct; % O1
P1.Bindex = 1;
P1.sPlocal = [-0.2; -0.25];

P2 = Point_struct; % O2
P2.Bindex = 2;
P2.sPlocal = [0.22; 0];

P3 = Point_struct; % A2
P3.Bindex = 2;
P3.sPlocal = [-0.22; 0];

P4 = Point_struct; % A3
P4.Bindex = 3;
P4.sPlocal = [-0.061501; -0.15310985];

P5 = Point_struct; % A2w
P5.Bindex = 2;
P5.sPlocal = [-0.22; 0];

P6 = Point_struct; % A5
P6.Bindex = 5;
P6.sPlocal = [0; 0];

P7 = Point_struct; % B3
P7.Bindex = 3;
P7.sPlocal = [0.061501; 0.15310985];

P8 = Point_struct; % B4
P8.Bindex = 4;
P8.sPlocal = [-0.17566585; 0.01189579];

P9 = Point_struct; % Q4
P9.Bindex = 4;
P9.sPlocal = [0.04133215; -0.02432392];

P10 = Point_struct; % Q1
P10.Bindex = 1;
P10.sPlocal = [-0.3; 0.02];

P11 = Point_struct; % C4
P11.Bindex = 4;
P11.sPlocal = [0.13433369; 0.01242813];

P12 = Point_struct; % D1
P12.Bindex = 1;
P12.sPlocal = [-0.18; -0.17];

P13 = Point_struct; % T1
P13.Bindex = 1;
P13.sPlocal = [0.32; 0.08];

P14 = Point_struct; % E6
P14.Bindex = 6;
P14.sPlocal = [0; 0];

Points = [P1; P2; P3; P4; P5; P6; P7; P8; P9; P10; P11; P12; P13; P14];
