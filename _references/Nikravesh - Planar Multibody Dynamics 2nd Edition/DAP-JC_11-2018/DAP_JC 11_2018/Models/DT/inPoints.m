function inPoints
    include_global

P1 = Point_struct; % A0
P1.Bindex = 0;
P1.sPlocal = [0; 0];

P2 = Point_struct; % B0
P2.Bindex = 0;
P2.sPlocal = [1.6; 0.2];

P3 = Point_struct; % A1
P3.Bindex = 1;
P3.sPlocal = [0.19696155; -0.03472964];

P4 = Point_struct; % C1
P4.Bindex = 1;
P4.sPlocal = [-0.19696155; 0.03472964];

P5 = Point_struct; % C2
P5.Bindex = 2;
P5.sPlocal = [-0.99534025; -0.096425];

P6 = Point_struct; % D2
P6.Bindex = 2;
P6.sPlocal = [0.29860208; 0.0289275];

P7 = Point_struct; % K2
P7.Bindex = 2;
P7.sPlocal = [0.99534025; 0.096425];

P8 = Point_struct; % D3
P8.Bindex = 3;
P8.sPlocal = [0.29999176; 0.00222352];

P9 = Point_struct; % E3
P9.Bindex = 3;
P9.sPlocal = [-0.99997253; -0.00741175];

P10 = Point_struct; % B3
P10.Bindex = 3;
P10.sPlocal = [0.99997253; 0.00741175];

P11 = Point_struct; % E4
P11.Bindex = 4;
P11.sPlocal = [-0.19847182; -0.02467662];

P12 = Point_struct; % F4
P12.Bindex = 4;
P12.sPlocal = [0.19847182; 0.02467662];

P13 = Point_struct; % F5
P13.Bindex = 5;
P13.sPlocal = [-0.8; -0.4];

P14 = Point_struct; % K5
P14.Bindex = 5;
P14.sPlocal = [0.8; -0.4];

Points = [P1; P2; P3; P4; P5; P6; P7; P8; P9; P10; P11; P12; P13; P14];
