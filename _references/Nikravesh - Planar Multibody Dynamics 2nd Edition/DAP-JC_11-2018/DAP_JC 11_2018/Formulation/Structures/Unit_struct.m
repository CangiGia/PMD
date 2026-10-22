function Unit = Unit_struct

Unit = struct ( ...
	'Bindex', 0	    ,   ... % body index
	'ulocal', [1;0] ,   ... % u_prime; xi and eta components
	'u'     , [0;0] ,	... % x, y components
	'u_r'   , [0;0] ,	... % vector u rotated
    'u_d'   , [0;0]     ... % u_dot
);
