function inPoints
    include_global

P1 = Point_struct; % O1
P1.Bindex = 1;
P1.sPlocal = [-0.25; -0.22];

P2 = Point_struct; % O2
P2.Bindex = 2;
P2.sPlocal = [0.22; 0];

P3 = Point_struct; % A2
P3.Bindex = 2;
P3.sPlocal = [-0.22; 0];

P4 = Point_struct; % A3
P4.Bindex = 3;
P4.sPlocal = [0; 0];

P5 = Point_struct; % B2
P5.Bindex = 2;
P5.sPlocal = [0.14090909; 0.19453696];

P6 = Point_struct; % C1
P6.Bindex = 1;
P6.sPlocal = [-0.15; -0.01];

P7 = Point_struct; % T1
P7.Bindex = 1;
P7.sPlocal = [0.26; 0.14];

P8 = Point_struct; % D4
P8.Bindex = 4;
P8.sPlocal = [0; 0];

Points = [P1; P2; P3; P4; P5; P6; P7; P8];
