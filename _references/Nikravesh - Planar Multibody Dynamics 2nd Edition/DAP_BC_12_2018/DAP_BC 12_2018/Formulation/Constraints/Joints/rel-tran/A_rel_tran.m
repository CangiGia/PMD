% rel-tran_A
% r-h-s of acc. constraint for relative-rotational constraint

    Pi = Joints(Ji).iPindex;  Pj = Joints(Ji).jPindex;
    Bi = Joints(Ji).iBindex;  Bj = Joints(Ji).jBindex;
    d  = Points(Pi).rP - Points(Pj).rP;
    d_d  = Points(Pi).rP_d - Points(Pj).rP_d;
    [fun, fun_d, fun_dd] = functs(Joints(Ji).iFunct, t);
    
    f = fun*fun_dd + fun_d^2;
    if Bi == 0
        f = f + d'*s_rot(Points(Pj).sP_d)'*Bodies(Bj).p_d;
    elseif Bj == 0
        f = f - d'*s_rot(Points(Pi).sP_d)'*Bodies(Bi).p_d - d_d'*d_d;
    else
        f = f + d'*s_rot(Points(Pj).sP_d)'*Bodies(Bj).p_d ...
              - d'*s_rot(Points(Pi).sP_d)'*Bodies(Bi).p_d - d_d'*d_d;
    end
