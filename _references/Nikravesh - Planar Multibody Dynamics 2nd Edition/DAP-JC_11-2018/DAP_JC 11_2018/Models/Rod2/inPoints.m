function inPoints
    include_global

P1 = Point_struct;
P1.Bindex = 1;
P1.sPlocal = [-0.5; 0];

P2 = Point_struct;
P2.Bindex = 1;
P2.sPlocal = [ 0.5; 0];

P3 = Point_struct;
P3.Bindex = 2;
P3.sPlocal = [-0.5; 0];

P4 = Point_struct;
P4.Bindex = 2;
P4.sPlocal = [ 0.5; 0];

Points = [P1; P2; P3; P4];
