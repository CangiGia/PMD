function inPoints
    include_global

A1 = Point_struct;
A1.Bindex = 1;
A1.sPlocal = [ 0.00; -0.07];

B1 = Point_struct;
B1.Bindex = 1;
B1.sPlocal = [-0.17;  0.25];

C1 = Point_struct;
C1.Bindex = 1;
C1.sPlocal = [ 0.11; -0.02];

O0 = Point_struct;
O0.Bindex = 0;
O0.sPlocal = [ 0.41;  0.83];

Q0 = Point_struct;
Q0.Bindex = 0;
Q0.sPlocal = [ 0.12; 0.29];

Points = [A1; B1; C1; O0; Q0];
