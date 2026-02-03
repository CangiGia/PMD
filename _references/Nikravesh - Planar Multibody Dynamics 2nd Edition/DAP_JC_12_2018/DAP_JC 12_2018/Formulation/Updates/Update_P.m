function Update_P(Pi)
% Recursively update coordinates of body i
    include_global
    
    Bi  = Points(Pi).Bindex;
    Points(Pi).sP = Bodies(Bi).A*Points(Pi).sPlocal;
    Points(Pi).sP_r = s_rot(Points(Pi).sP);
    Points(Pi).rP = Bodies(Bi).r + Points(Pi).sP;
    