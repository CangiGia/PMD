% rev_tran_C
% Evaluates the constraint for a revolute-translational joint 

    Pi = Joints(Ji).iPindex;  Pj = Joints(Ji).jPindex;
    ui_r = Uvectors(Joints(Ji).iUindex).u_r;
    d  = Points(Pi).rP - Points(Pj).rP;
        f = ui_r'*d - Joints(Ji).L;
