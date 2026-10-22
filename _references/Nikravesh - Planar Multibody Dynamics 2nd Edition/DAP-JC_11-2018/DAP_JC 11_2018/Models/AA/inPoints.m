function inPoints
    include_global

Q1 = Point_struct;
Q1.Bindex = 1;
Q1.sPlocal = [-0.24; 0];

A1 = Point_struct;
A1.Bindex = 1;
A1.sPlocal = [0.18;0];

A2 = Point_struct;
A2.Bindex = 2;
A2.sPlocal = [-0.07;-0.10];

B2 = Point_struct;
B2.Bindex = 2;
B2.sPlocal = [-0.10;0.12];

B3 = Point_struct;
B3.Bindex = 3;
B3.sPlocal = [0.13;0];

O3 = Point_struct;
O3.Bindex = 3;
O3.sPlocal = [-0.13;0];

O0 = Point_struct;
O0.Bindex = 0;
O0.sPlocal = [0.32;0.40];

Q0 = Point_struct;
Q0.Bindex = 0;
Q0.sPlocal = [0.20;0.26];

E1 = Point_struct;
E1.Bindex = 1;
E1.sPlocal = [0;0];

F0 = Point_struct;
F0.Bindex = 0;
F0.sPlocal = [0.38;0.43];

Points = [Q1; A1; A2; B2; B3; O3; O0; Q0; E1; F0];
