function Joint = Joint_struct

Joint = struct ( ...
	'type'     , 'rev' , ... % joint type: rev, tran, rev-rev, rev-tran, rigid, disc, rel-rot, rel-tran
    'iBindex'  , 0     , ... % body index i
    'jBindex'  , 0     , ... % body index j
    'iPindex'  , 0     , ... % point Pi index
    'jPindex'  , 0     , ... % point Pj index
    'iUindex'  , 0     , ... % unit vector u_i index
    'jUindex'  , 0     , ... % unit vector u_j index
    'iFunct'   , 0     , ... % analytical function index
	'L' 	   , 0     , ... % constant length
	'R' 	   , 1     , ... % constant radius
    'x0'       , 0     , ... % initial condition x for disc
    'p0'       , 0     , ... % initial condition phi for a disc (or rigid)
    'd0'       , []    , ... % initial condition for d (rigid)
    'fix'      , 0     , ... % fix relative dof if = 1 (rev or tran)
    'nbody'    , 2     , ... % number of moving bodies involved
    'mrows'    , 2     , ... % number of rows (constraints)
    'rows'     , 0     , ... % row index-start
    'rowe'     , 0     , ... % row index-end
    'colis'    , 0     , ... % column index for body i-start
    'colie'    , 0     , ... % column index for body i-end
    'coljs'    , 0     , ... % column index for body j-start
    'colje'    , 0     , ... % column index for body j-end
    'lagrange' , zeros(3,1) ... % Lagrange multipliers
);
