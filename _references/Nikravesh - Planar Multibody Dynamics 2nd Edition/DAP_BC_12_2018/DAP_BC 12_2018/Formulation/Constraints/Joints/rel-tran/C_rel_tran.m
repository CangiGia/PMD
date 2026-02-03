% rel-tran_C
% Constraint for relative-translational 

    Pi = Joints(Ji).iPindex;  Pj = Joints(Ji).jPindex;
    d  = Points(Pi).rP - Points(Pj).rP;
    [fun, fun_d, fun_dd] = functs(Joints(Ji).iFunct, t);
    
        f = (d'*d - fun^2)/2;
