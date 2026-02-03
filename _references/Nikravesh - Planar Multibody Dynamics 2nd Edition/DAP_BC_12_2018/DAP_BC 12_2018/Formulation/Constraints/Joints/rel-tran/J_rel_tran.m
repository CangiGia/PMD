% rel-tran_J
% Jacobian sub-matrices for relative-translational constraint 

    Pi = Joints(Ji).iPindex;  Pj = Joints(Ji).jPindex;
    d  = Points(Pi).rP - Points(Pj).rP;
        Di = [ d'  d'*Points(Pi).sP_r];
        Dj = [-d' -d'*Points(Pj).sP_r];
