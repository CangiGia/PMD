function user_force
    include_global

% Parameters for the contact model
    k = 10^11; e = 0.95;
    
% Contact(Contact index, Point index, Body index, k, e, Model index)    
%    Model-index = 1: Eq. 11.42
%    Model-index = 2: Eq. 11.43

% Point [1] on Body (1)
    Contact(1, 1, 1, k, e, 1)
% Point [2] on Body (1)
    Contact(2, 2, 1, k, e, 1)
 % Point [4] on Body (2)
    Contact(3, 4, 2, k, e, 1)
    
