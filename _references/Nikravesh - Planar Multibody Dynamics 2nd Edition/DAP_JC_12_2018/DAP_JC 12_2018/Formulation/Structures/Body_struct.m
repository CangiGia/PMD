function Body = Body_struct

Body = struct ( ...
    'm'     , 1     , ... % mass
    'J'     , 1     , ... % moment of inertia 
	'r'	    , [0;0]	, ... % x, y coordinates
	'p'  	, 0     , ... % angle phi
	'r_d'   , [0;0]	, ... % time derivative of x and y
    'p_d'   , 0     , ... % time derivative of phi
	'A'	    , eye(2), ... % rotational transformation matrix
	'r_dd'  , [0;0]	, ... % x_double_dot,y_double_do
	'p_dd'  , 0   	, ... % 2nd time derivative of phi
    'theta' , 0     , ... % joint coordinate(s)
    'theta_d', 0    , ... % joint velocity(ies) 
    'nJC'   , 1     , ... % number of joint coordinates
    'p0'    , 0     , ... % constant angle for tran joint (phi_i -phi_im)
    'irc'   , 0     , ... % index of the 1st element of r in u or u_dot
    'm_inv' , 1     , ... % mass inverse
    'J_inv' , 1     , ... % inverse of moment of inertia
    'wgt'   , [0;0] , ... % weight of body as a force vector
    'f'     , [0;0] , ... % sum of forces that act on the body
    'n'     , 0     , ... % sum of moments that act on the body
    'shape' , ' '   , ... % 'circle', 'rect', line
    'R'     , 1     , ... % radius of the circle
    'circ'  , []    , ... % points on circumference of the circle
    'W'     , 0     , ... % width of the rectangle
    'H'     , 0     , ... % hight of the rectangle
    'color' , 'k'   , ... % default color for the body
    'P4'    , []    , ... % 4 corners of the rectangle
    'pts'   , []      ... % point indexes associated with this body
);

