function Coord_Tran(Pim, Uim, Pi, Ui, theta)
% Recursively update coordinates of body i
    include_global
    
    Bim = Points(Pim).Bindex;
    Bi  = Points(Pi).Bindex;
    Bodies(Bi).theta = theta;
    if Bim == 0
        Bodies(Bi).p = Bodies(Bi).p0;
    else
        Bodies(Bi).p = Bodies(Bim).p + Bodies(Bi).p0;
        Points(Pim).sP = Bodies(Bim).A * Points(Pim).sPlocal;
        Points(Pim).sP_r = s_rot(Points(Pim).sP);
        Uvectors(Uim).u = Bodies(Bim).A * Uvectors(Uim).ulocal;
        Uvectors(Uim).u_r = s_rot(Uvectors(Uim).u);
    end
    Points(Pi).rP = Points(Pim).rP + theta*Uvectors(Uim).u;
    Uvectors(Ui).u = Uvectors(Uim).u;
    Uvectors(Ui).u_r = s_rot(Uvectors(Ui).u);    
    Bodies(Bi).A    = Matrix_A(Bodies(Bi).p);
    Points(Pi).sP   = Bodies(Bi).A*Points(Pi).sPlocal;
    Points(Pi).sP_r = s_rot(Points(Pi).sP);
    Bodies(Bi).r    = Points(Pi).rP - Points(Pi).sP;