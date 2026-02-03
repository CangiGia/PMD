function inPoints
    include_global

P1 = Point_struct;
P1.Bindex = 0;
P1.sPlocal = [ 0; 0];

P2 = Point_struct;
P2.Bindex = 1;
P2.sPlocal = [ 0; 0.2];

P3 = Point_struct;
P3.Bindex = 1;
P3.sPlocal = [-0.15; 0.2];

P4 = Point_struct;
P4.Bindex = 1;
P4.sPlocal = [ 0.15; 0.2];

P5 = Point_struct;
P5.Bindex = 2;
P5.sPlocal = [ 0; 0.15];

P6 = Point_struct;
P6.Bindex = 3;
P6.sPlocal = [ 0; 0.15];

Points = [P1; P2; P3; P4; P5; P6];
