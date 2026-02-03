% rel-rot_C
% Constraint for relative-rotational 

    [fun, fun_d, fun_dd] = functs(Joints(Ji).iFunct, t);
    Bi = Joints(Ji).iBindex; Bj = Joints(Ji).jBindex;
    
    if Bi == 0
        f = -Bodies(Bj).p - fun;
    elseif Bj == 0
        f =  Bodies(Bi).p - fun;
    else
        f =  Bodies(Bi).p - Bodies(Bj).p - fun;
    end
