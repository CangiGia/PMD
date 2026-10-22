function Update_U(Ui)
% Recursively update coordinates of body i
    include_global
    
    Bi  = Uvectors(Ui).Bindex;
    Uvectors(Ui).u = Bodies(Bi).A*Uvectors(Ui).ulocal;
    Uvectors(Ui).u_r = s_rot(Uvectors(Ui).u);
    