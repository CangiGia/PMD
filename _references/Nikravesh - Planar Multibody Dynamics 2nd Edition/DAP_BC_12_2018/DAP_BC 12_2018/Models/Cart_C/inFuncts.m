function inFuncts
    include_global

F1 = Funct_struct;
F1.type = 'a';
F1.coeff = [0 -2*pi 0];

% F1.type = 'c';
% F1.t_end = 2.0;
% F1.dfdt_end = -2*pi;

Functs = [F1];
