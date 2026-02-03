% revolute_translational_A
% r-h-s of acc. constraint for a revolute translational joint

    Pi = Joints(Ji).iPindex;  Pj = Joints(Ji).jPindex;
    Bi = Joints(Ji).iBindex;  Bj = Joints(Ji).jBindex;
    ui  = Uvectors(Joints(Ji).iUindex).u;
    ui_d = Uvectors(Joints(Ji).iUindex).u_d;
    d  = Points(Pi).rP - Points(Pj).rP;
    d_d  = Points(Pi).rP_d - Points(Pj).rP_d;
    
    if Bi == 0
        f = ui'*Points(Pj).sP_d*Bodies(Bj).p_d;
    elseif Bj == 0
        f = ui_d'*(d*Bodies(Bi).p_d + 2*s_rot(d_d)) - ...
            ui'*Points(Pi).sP_d*Bodies(Bi).p_d;
    else
        f = ui_d'*(d*Bodies(Bi).p_d + 2*s_rot(d_d)) - ...
            ui'*(Points(Pi).sP_d*Bodies(Bi).p_d - ...
                  Points(Pj).sP_d*Bodies(Bj).p_d);
    end
