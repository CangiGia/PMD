function Funct = Funct_struct

Funct = struct ( ...
	'type'     , 'a'   , ... % function type a, b, or c
    't_start'  , 0     , ... % required for functions b, c
    'f_start'  , 0     , ... % required for functions b, c
    't_end'    , 1     , ... % required for functions b, c
    'f_end'    , 1     , ... % required for functions b
    'dfdt_end' , 1     , ... % required for functions c
    'ncoeff'   , 4     , ... % number of coefficients
    'coeff'    , []      ... % required for function a
);
