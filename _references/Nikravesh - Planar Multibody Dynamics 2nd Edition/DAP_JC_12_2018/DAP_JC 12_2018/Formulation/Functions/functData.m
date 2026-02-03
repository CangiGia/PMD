 function functData(Ci)
    include_global

switch (Functs(Ci).type)
    case {'a'}
        Functs(Ci).ncoeff = 4;
        Functs(Ci).coeff(4) = 2*Functs(Ci).coeff(3);
    case {'b'}
        Functs(Ci).ncoeff = 9;
        xe = Functs(Ci).t_end - Functs(Ci).t_start;
        fe = Functs(Ci).f_end - Functs(Ci).f_start;
        C = [  xe^3     xe^4      xe^5
             3*xe^2   4*xe^3    5*xe^4
             6*xe    12*xe^2   20*xe^3];
        sol = C\[fe; 0; 0];
        Functs(Ci).coeff(1:3) = sol';
        Functs(Ci).coeff(4) =  3*sol(1);
        Functs(Ci).coeff(5) =  4*sol(2);
        Functs(Ci).coeff(6) =  5*sol(3);
        Functs(Ci).coeff(7) =  6*sol(1);
        Functs(Ci).coeff(8) = 12*sol(2);
        Functs(Ci).coeff(9) = 20*sol(3);
    case {'c'}
        Functs(Ci).ncoeff = 9;
        xe = Functs(Ci).t_end - Functs(Ci).t_start;
        fpe = Functs(Ci).dfdt_end;
        C = [4*xe^3   5*xe^4    6*xe^5
            12*xe^2  20*xe^3   30*xe^4
            24*xe    60*xe^2  120*xe^3];
        sol = C\[fpe; 0; 0];
        Functs(Ci).coeff(1:3) = sol';
        Functs(Ci).coeff(4) =  4*sol(1);
        Functs(Ci).coeff(5) =  5*sol(2);
        Functs(Ci).coeff(6) =  6*sol(3);
        Functs(Ci).coeff(7) = 12*sol(1);
        Functs(Ci).coeff(8) = 20*sol(2);
        Functs(Ci).coeff(9) = 30*sol(3);
end


