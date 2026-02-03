function inPoints
    include_global

P1 = Point_struct;
P1.Bindex = 0;
P1.sPlocal = [ 0; 0.2];

P2 = Point_struct;
P2.Bindex = 1;
P2.sPlocal = [ 0; 0];

Points = [P1; P2];
