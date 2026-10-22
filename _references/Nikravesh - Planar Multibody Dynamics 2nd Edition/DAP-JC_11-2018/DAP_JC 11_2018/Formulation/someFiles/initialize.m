function initialize
    include_global
    global theta theta_d
    bodycolor = {'r' 'g' 'b' 'c' 'm'};

    num = 0; % number of function evaluations
    t10 = 0;
    flags = zeros(10,1); pen_d0 = zeros(10,1);
    
% Bodies
    nB = length(Bodies); nB3 = 3*nB; 
    nJC = 0; theta = []; theta_d = [];
    for Bi = 1:nB
        Bodies(Bi).irc    = 3*(Bi-1) + 1;
        Bodies(Bi).m_inv  = 1/Bodies(Bi).m;
        Bodies(Bi).J_inv  = 1/Bodies(Bi).J;
        Bodies(Bi).nJCB   = length(Bodies(Bi).theta);
        if Bodies(Bi).nJCB > 1
            Bodies(Bi).theta = Bodies(Bi).theta(:);
            Bodies(Bi).theta_d = Bodies(Bi).theta_d(:);
        end
        nJC = nJC + Bodies(Bi).nJCB;
        nTemp = length(Bodies(Bi).theta_d);
        if nTemp ~= Bodies(Bi).nJCB
            Bodies(Bi).theta_d = [0 0 0]';
        end
        theta   = [theta;   Bodies(Bi).theta];
        theta_d = [theta_d; Bodies(Bi).theta_d];
    end
% Theta's
    nJC2 = 2*nJC;    
% Mass (inertia) matrix as an array
    M_array = zeros(nB3,1); M_inv_array = zeros(nB3,1);
    for Bi = 1:nB
        is = 3*(Bi - 1) + 1;
        ie = is + 2;
        M_array(is:ie,1) = [Bodies(Bi).m; Bodies(Bi).m; Bodies(Bi).J];
        M_inv_array(is:ie,1) = [Bodies(Bi).m_inv; Bodies(Bi).m_inv; ...
                                Bodies(Bi).J_inv];
    end
    
% Points
    nP = length(Points); nPanim = length(Points_anim);
    nPtot = nP + nPanim;
    Points = [Points; Points_anim];
    for Pi = 1:nPtot
        if Points(Pi).Bindex == 0
            Points(Pi).sP   = Points(Pi).sPlocal;
            Points(Pi).sP_r = s_rot(Points(Pi).sP);
            Points(Pi).rP   = Points(Pi).sP;
        end
        for Bi = 1:nB
            if Points(Pi).Bindex == Bi
                len = length(Bodies(Bi).pts); %current length of pts
                Bodies(Bi).pts(len + 1) = Pi;
            end
        end
    end
    
% Unit vectors  
   nU = length(Uvectors); 
    for Vi = 1:nU
        if Uvectors(Vi).Bindex == 0
            Uvectors(Vi).u = Uvectors(Vi).ulocal;
            Uvectors(Vi).u_r = s_rot(Uvectors(Vi).u);
        end
    end

% Force elements
    nF = length(Forces);
    for Fi = 1:nF
        switch (Forces(Fi).type);
            case {'weight'}
                ug = Forces(Fi).gravity*Forces(Fi).wgt;
                for Bi = 1:nB
                    Bodies(Bi).wgt = Bodies(Bi).m*ug;
                end
            case {'ptp'}
                Pi = Forces(Fi).iPindex;    Pj = Forces(Fi).jPindex;
                Forces(Fi).iBindex = Points(Pi).Bindex;
                Forces(Fi).jBindex = Points(Pj).Bindex;
        end
    end
                        
% Functions
    nFc = length(Functs);
    for Ci = 1:nFc
        functData(Ci);
    end

