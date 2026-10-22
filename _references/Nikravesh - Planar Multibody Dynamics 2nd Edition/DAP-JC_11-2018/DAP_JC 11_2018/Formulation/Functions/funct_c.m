function [f, f_d, f_dd] = funct_c(Ci, xx)
% Function type 'c'
    include_global

c = Functs(Ci).coeff;
if xx <= Functs(Ci).t_start
    f    = Functs(Ci).f_start;
    f_d  = 0;
    f_dd = 0;
elseif xx > Functs(Ci).t_start && xx < Functs(Ci).t_end
    x    = xx - Functs(Ci).t_start;
    f    = c(1)*x^4 + c(2)*x^5 + c(3)*x^6 + Functs(Ci).f_start;
    f_d  = c(4)*x^3 + c(5)*x^4 + c(6)*x^5;
    f_dd = c(7)*x^2 + c(8)*x^3 + c(9)*x^4;
else
    f    = 0; % this should be undefined
    f_d  = Functs(Ci).dfdt_end;
    f_dd = 0;
end
