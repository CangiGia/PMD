% translational_J
% Jacobian sub-matrices for a translational joint

    Pi = Joints(Ji).iPindex;      Pj = Joints(Ji).jPindex;
    uj = Uvectors(Joints(Ji).jUindex).u; 
    uj_r = Uvectors(Joints(Ji).jUindex).u_r;
    d  = Points(Pi).rP - Points(Pj).rP;
    
        Di = [ uj_r'  uj'*Points(Pi).sP
               0 0       1];
        Dj = [-uj_r' -uj'*(Points(Pj).sP + d)
               0 0      -1];
        
    if Joints(Ji).fix == 1
        Di = [Di
              uj'  uj'*Points(Pi).sP_r];
        Dj = [Dj
              -uj' -uj'*Points(Pj).sP_r];
    end           
 