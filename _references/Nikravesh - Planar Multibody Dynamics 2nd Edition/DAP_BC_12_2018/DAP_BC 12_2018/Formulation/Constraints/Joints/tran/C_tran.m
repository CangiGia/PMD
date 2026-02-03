% translational_C
% Evaluate the constraints for a translational joint 

    Pi = Joints(Ji).iPindex;    Pj = Joints(Ji).jPindex;

    uj_r = Uvectors(Joints(Ji).jUindex).u_r;
    ui = Uvectors(Joints(Ji).iUindex).u;
    d  = Points(Pi).rP - Points(Pj).rP;
    
    f  = [uj_r'*d; uj_r'*ui];
    
    if Joints(Ji).fix == 1
        f = [f
            (ui'*d - Joints(Ji).p0)/2];
    end

