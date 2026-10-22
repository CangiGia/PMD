function inPoints
    include_global

P1 = Point_struct; % G0
P1.Bindex = 0;
P1.sPlocal = [0; 0.05];

P2 = Point_struct; % Tc
P2.Bindex = 1;
P2.sPlocal = [0; 0];

P3 = Point_struct; % R0a
P3.Bindex = 1;
P3.sPlocal = [0; 0.03];

P4 = Point_struct; % R0b
P4.Bindex = 2;
P4.sPlocal = [-0.00236474; -0.0075];

P5 = Point_struct; % R1a
P5.Bindex = 2;
P5.sPlocal = [0.00236474; 0.0075];

P6 = Point_struct; % R1b
P6.Bindex = 3;
P6.sPlocal = [-0.00236474; -0.0075];

P7 = Point_struct; % R2a
P7.Bindex = 3;
P7.sPlocal = [0.00236474; 0.0075];

P8 = Point_struct; % R2b
P8.Bindex = 4;
P8.sPlocal = [-0.00236474; -0.0075];

P9 = Point_struct; % R3a
P9.Bindex = 4;
P9.sPlocal = [0.00236474; 0.0075];

P10 = Point_struct; % R3b
P10.Bindex = 5;
P10.sPlocal = [-0.00236474; -0.0075];

P11 = Point_struct; % R4a
P11.Bindex = 5;
P11.sPlocal = [0.00236474; 0.0075];

P12 = Point_struct; % R4b
P12.Bindex = 6;
P12.sPlocal = [-0.00236474; -0.0075];

P13 = Point_struct; % R5a
P13.Bindex = 6;
P13.sPlocal = [0.00236474; 0.0075];

P14 = Point_struct; % R5b
P14.Bindex = 7;
P14.sPlocal = [-0.00236474; -0.0075];

P15 = Point_struct; % R6a
P15.Bindex = 7;
P15.sPlocal = [0.00236474; 0.0075];

P16 = Point_struct; % R6b
P16.Bindex = 8;
P16.sPlocal = [-0.00236474; -0.0075];

P17 = Point_struct; % R7a
P17.Bindex = 8;
P17.sPlocal = [0.00236474; 0.0075];

P18 = Point_struct; % R7b
P18.Bindex = 9;
P18.sPlocal = [-0.033; -0.046];

P19 = Point_struct; % Sk1
P19.Bindex = 9;
P19.sPlocal = [0.02; 0.02];

P20 = Point_struct; % Sk2
P20.Bindex = 9;
P20.sPlocal = [-0.06; 0.02];

P21 = Point_struct; % Sk3
P21.Bindex = 9;
P21.sPlocal = [-0.02; 0.06];

P22 = Point_struct; % Sk4
P22.Bindex = 9;
P22.sPlocal = [-0.02; -0.02];

P23 = Point_struct; % Bc
P23.Bindex = 10;
P23.sPlocal = [0; 0];

Points = [P1; P2; P3; P4; P5; P6; P7; P8; P9; P10; P11; P12; P13; P14; P15; P16; P17; P18; P19; P20; P21; P22; P23];
