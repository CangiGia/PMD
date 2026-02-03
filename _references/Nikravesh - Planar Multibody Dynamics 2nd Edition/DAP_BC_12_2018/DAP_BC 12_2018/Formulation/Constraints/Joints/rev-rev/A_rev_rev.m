% revolute_revolute_A
% r-h-s of acc. constraint for a revolute_revolute joint

    Pi = Joints(Ji).iPindex;  Pj = Joints(Ji).jPindex;
    Bi = Joints(Ji).iBindex;  Bj = Joints(Ji).jBindex;
    d  = Points(Pi).rP - Points(Pj).rP;
    d_d  = Points(Pi).rP_d - Points(Pj).rP_d;

    L = Joints(Ji).L; u = d/L; u_d = d_d/L;
    
        f = - u_d'*d_d;
    if Bi == 0
        f = f + u'*s_rot(Points(Pj).sP_d)*Bodies(Bj).p_d;
    elseif Bj == 0
        f = f - u'*s_rot(Points(Pi).sP_d)*Bodies(Bi).p_d;
    else
        f = f - u'*(s_rot(Points(Pi).sP_d*Bodies(Bi).p_d - ...
                          Points(Pj).sP_d*Bodies(Bj).p_d));
    end
