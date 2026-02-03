function [f, f_d, f_dd] = funct_a(Ci, x)
% Function type 'a'
    include_global

    c = Functs(Ci).coeff;

    f    = c(1) + c(2)*x + c(3)*x^2;
    f_d  = c(2) + c(4)*x;
    f_dd = c(4);
