function Coord_Float(Bi, theta)
% Recursively update coordinates of body i
    include_global
    
    Bodies(Bi).r = theta(1:2);
    Bodies(Bi).p = theta(3);
    Bodies(Bi).theta = theta;
    Bodies(Bi).A = Matrix_A(Bodies(Bi).p);
    