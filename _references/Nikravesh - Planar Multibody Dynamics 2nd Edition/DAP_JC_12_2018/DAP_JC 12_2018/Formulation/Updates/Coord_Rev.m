function Coord_Rev(Pim, Pi, theta)
% Recursively update coordinates of body i
    include_global
    
    Bim = Points(Pim).Bindex;
    Bi  = Points(Pi).Bindex;
    Bodies(Bi).theta = theta;
    if Bim == 0
        Bodies(Bi).p = Bodies(Bi).theta;
    else
        Points(Pim).sP   = Bodies(Bim).A*Points(Pim).sPlocal;
        Points(Pim).sP_r = s_rot(Points(Pim).sP);
        Points(Pim).rP   = Bodies(Bim).r + Points(Pim).sP;
        Bodies(Bi).p     = Bodies(Bim).p + theta;
    end
    Points(Pi).rP   = Points(Pim).rP;
    Bodies(Bi).A    = Matrix_A(Bodies(Bi).p);
    Points(Pi).sP   = Bodies(Bi).A*Points(Pi).sPlocal;
    Points(Pi).sP_r = s_rot(Points(Pi).sP);
    Bodies(Bi).r    = Points(Pi).rP - Points(Pi).sP;
