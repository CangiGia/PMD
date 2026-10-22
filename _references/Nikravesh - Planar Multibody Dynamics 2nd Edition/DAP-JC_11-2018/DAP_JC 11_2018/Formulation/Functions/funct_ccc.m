function [f, f_d, f_dd] = funct_c(Ci, xx)
% Function type 'c'
    include_global

c = Functs(Ci).coeff;
if xx <= Functs(Ci).x_start
    f = Functs(Ci).f_start;
    f_d = 0;
    f_dd = 0;
elseif xx > Functs(Ci).x_start && xx < Functs(Ci).x_end
    x = xx - Functs(Ci).x_start;
    f    =  c(1)*x^3 +  c(2)*x^4 +  c(3)*x^5 +  c(4)*x^6 +  c(5)*x^7 + Functs(Ci).f_start;
    f_d  =  c(6)*x^2 +  c(7)*x^3 +  c(8)*x^4 +  c(9)*x^5 + c(10)*x^6;
    f_dd = c(11)*x   + c(12)*x^2 + c(13)*x^3 + c(14)*x^4 + c(15)*x^5;
else
    f = Functs(Ci).f_end;
    f_d = 0;
    f_dd = 0;
end
