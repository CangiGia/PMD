% rel-rot_V
% r-h-s of vel. constraint for relative-rotational constraint

    [fun, fun_d, fun_dd] = functs(Joints(Ji).iFunct, t);
    f = fun_d;
