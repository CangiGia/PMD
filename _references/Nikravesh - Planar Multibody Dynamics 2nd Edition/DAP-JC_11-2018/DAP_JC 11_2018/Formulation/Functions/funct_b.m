function [f, f_d, f_dd] = funct_b(Ci, xx)
% Function type 'b'
    include_global

c = Functs(Ci).coeff;
if xx <= Functs(Ci).t_start
    f    = Functs(Ci).f_start;
    f_d  = 0;
    f_dd = 0;
elseif xx > Functs(Ci).t_start && xx < Functs(Ci).t_end
    x    = xx - Functs(Ci).t_start;
    f    = c(1)*x^3 + c(2)*x^4 + c(3)*x^5 + Functs(Ci).f_start;
    f_d  = c(4)*x^2 + c(5)*x^3 + c(6)*x^4;
    f_dd = c(7)*x   + c(8)*x^2 + c(9)*x^3;
else
    f    = Functs(Ci).f_end;
    f_d  = 0;
    f_dd = 0;
end
