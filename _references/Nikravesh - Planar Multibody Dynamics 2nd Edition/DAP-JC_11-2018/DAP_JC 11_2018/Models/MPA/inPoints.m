function inPoints
    include_global
    
Q1 = Point_struct;
Q1.Bindex = 1;
Q1.sPlocal = [-0.225; 0.00];

A1 = Point_struct;
A1.Bindex = 1;
A1.sPlocal = [ 0.225; 0.00];

A2 = Point_struct;
A2.Bindex = 2;
A2.sPlocal = [ 0.00; -0.07];

B2 = Point_struct;
B2.Bindex = 2;
B2.sPlocal = [-0.17;  0.25];

C2 = Point_struct;
C2.Bindex = 2;
C2.sPlocal = [ 0.11; -0.02];

O3 = Point_struct;
O3.Bindex = 3;
O3.sPlocal = [-0.15;  0.00];

O0 = Point_struct;
O0.Bindex = 0;
O0.sPlocal = [ 0.41;  0.83];

Q0 = Point_struct;
Q0.Bindex = 0;
Q0.sPlocal = [ 0.12; 0.29];

Points = [Q1; A1; A2; B2; C2; O3; O0; Q0];
