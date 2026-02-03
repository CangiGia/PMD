% rev_tran_J
% Jacobian sub-matrices for a revolute-translational joint

    Pi = Joints(Ji).iPindex;  Pj = Joints(Ji).jPindex;
    ui = Uvectors(Joints(Ji).iUindex).u; 
    ui_r = Uvectors(Joints(Ji).iUindex).u_r;
    d  = Points(Pi).rP - Points(Pj).rP;
        Di = [ ui_r'  ui'*(Points(Pi).sP - d)];
        Dj = [-ui_r' -ui'*Points(Pj).sP];