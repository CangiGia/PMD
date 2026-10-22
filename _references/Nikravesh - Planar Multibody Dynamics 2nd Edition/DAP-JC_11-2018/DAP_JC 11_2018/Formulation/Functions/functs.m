function [f, f_d, f_dd] = functs(Ci, t)
    include_global

    switch (Functs(Ci).type)
        case {'a'}
            [f, f_d, f_dd] = funct_a(Ci, t);
        case {'b'}
            [f, f_d, f_dd] = funct_b(Ci, t);
        case {'c'}
            [f, f_d, f_dd] = funct_c(Ci, t);
        case {'d'}
            [f, f_d, f_dd] = funct_d(Ci, t);
    end