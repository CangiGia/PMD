function Point = Point_struct

Point = struct ( ...
	'Bindex' , 0     , ... % body index
	'sPlocal', [0;0] , ... % body-fixed coordinates
	'sP'     , [0;0] , ... % x, y components of vector s
	'sP_r'   , [0;0] , ... % vector s rotated
	'rP'     , [0;0] , ... % x, y coordinates of the point
    'sP_d'   , [0;0] , ... % s_P_dot
    'rP_d'   , [0;0] , ... % r_P_dot
    'rP_dd'  , [0;0]   ... % r_P_dot2
);
