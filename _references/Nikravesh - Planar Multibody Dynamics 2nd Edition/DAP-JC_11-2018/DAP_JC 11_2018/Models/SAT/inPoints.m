function inPoints
    include_global

P1 = Point_struct; % B0
P1.Bindex = 0;
P1.sPlocal = [-0.6; -0.04];

P2 = Point_struct; % C0
P2.Bindex = 0;
P2.sPlocal = [-0.3; -0.02];

P3 = Point_struct; % D0
P3.Bindex = 0;
P3.sPlocal = [0; -0.2];

P4 = Point_struct; % B1
P4.Bindex = 1;
P4.sPlocal = [0.3; -0.04];

P5 = Point_struct; % A1
P5.Bindex = 1;
P5.sPlocal = [0; -0.02];

P6 = Point_struct; % A2
P6.Bindex = 2;
P6.sPlocal = [-0.175; 0];

P7 = Point_struct; % S12
P7.Bindex = 2;
P7.sPlocal = [0.175; 0];

P8 = Point_struct; % S13
P8.Bindex = 3;
P8.sPlocal = [-0.175; 0];

P9 = Point_struct; % S13b
P9.Bindex = 3;
P9.sPlocal = [-0.175; 0];

P10 = Point_struct; % C3
P10.Bindex = 3;
P10.sPlocal = [0.175; 0];

P11 = Point_struct; % S14
P11.Bindex = 4;
P11.sPlocal = [-0.15; 0];

P12 = Point_struct; % L14
P12.Bindex = 4;
P12.sPlocal = [0.15; 0];

P13 = Point_struct; % L15
P13.Bindex = 5;
P13.sPlocal = [-0.15; 0];

P14 = Point_struct; % D5
P14.Bindex = 5;
P14.sPlocal = [0.15; 0];

Points = [P1; P2; P3; P4; P5; P6; P7; P8; P9; P10; P11; P12; P13; P14];
