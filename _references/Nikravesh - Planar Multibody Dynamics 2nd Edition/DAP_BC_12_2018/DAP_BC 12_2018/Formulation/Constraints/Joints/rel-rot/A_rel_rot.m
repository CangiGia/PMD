% rel-rot_A
% r-h-s of acc. constraint for relative-rotational constraint

    [fun, fun_d, fun_dd] = functs(Joints(Ji).iFunct, t);
    f = fun_dd;
