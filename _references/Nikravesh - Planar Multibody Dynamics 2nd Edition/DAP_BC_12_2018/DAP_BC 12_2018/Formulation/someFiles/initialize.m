function initialize
include_global
bodycolor = {'r' 'g' 'b' 'c' 'm'};

    num = 0; % number of function evaluations
    t10 = 0;
    flags = zeros(10,1); pen_d0 = zeros(10,1);
    
% Bodies
    nB = length(Bodies); nB3 = 3*nB; nB6 = 6*nB;
    for Bi = 1:nB
        Bodies(Bi).irc    = 3*(Bi-1) + 1;
        Bodies(Bi).irv    = nB3 + 3*(Bi-1) + 1;
        Bodies(Bi).m_inv  = 1/Bodies(Bi).m;
        Bodies(Bi).J_inv  = 1/Bodies(Bi).J;
        Bodies(Bi).A      = Matrix_A(Bodies(Bi).p);
%         Bodies(Bi).color  = bodycolor{mod(Bi-1, 5) + 1};
    end
    
% Mass (inertia) matrix as an array
    M_array = zeros(nB3,1); M_inv_array = zeros(nB3,1);
    for Bi = 1:nB
        is = 3*(Bi - 1) + 1;
        ie = is + 2;
        M_array(is:ie,1) = [Bodies(Bi).m; Bodies(Bi).m; Bodies(Bi).J];
        M_inv_array(is:ie,1) = [Bodies(Bi).m_inv; Bodies(Bi).m_inv; Bodies(Bi).J_inv];
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
                        
% Joints
   nJ = length(Joints); 
   cfriction = 0;
    % Assign number of constraints and number of bodies to each joint type
    
for Ji = 1:nJ
    switch (Joints(Ji).type);
        case {'rev'};
            Joints(Ji).mrows = 2; Joints(Ji).nbody = 2;
            Pi = Joints(Ji).iPindex;    Pj = Joints(Ji).jPindex;
            Bi = Points(Pi).Bindex; Joints(Ji).iBindex = Bi;
            Bj = Points(Pj).Bindex; Joints(Ji).jBindex = Bj;
            if Joints(Ji).fix == 1
                Joints(Ji).mrows = 3;
%                 Joints(Ji).p0 = Bodies(Bi).p - Bodies(Bj).p;
                if Bi == 0
                    Joints(Ji).p0 = - Bodies(Bj).p;
                elseif Bj == 0
                    Joints(Ji).p0 = Bodies(Bi).p;
                else
                    Joints(Ji).p0 = Bodies(Bi).p - Bodies(Bj).p;
                end
            end
        case {'tran'}
            Joints(Ji).mrows = 2; Joints(Ji).nbody = 2;
            Pi = Joints(Ji).iPindex; Pj = Joints(Ji).jPindex;
            Bi = Points(Pi).Bindex; Joints(Ji).iBindex = Bi;
            Bj = Points(Pj).Bindex; Joints(Ji).jBindex = Bj;
            if Joints(Ji).fix == 1
                Joints(Ji).mrows = 3;
                if Bi == 0
                    Joints(Ji).p0 = norm(Points(Pi).rP - ...
                        Bodies(Bj).r - Bodies(Bj).A*Points(Pj).sPlocal);
                elseif Bj == 0
                    Joints(Ji).p0 = norm(Bodies(Bi).r + ...
                        Bodies(Bi).A*Points(Pi).sPlocal - Points(Pj).rP);
                else
                    Joints(Ji).p0 = norm(Bodies(Bi).r + ...
                        Bodies(Bi).A*Points(Pi).sPlocal - ...
                        Bodies(Bj).r - Bodies(Bj).A*Points(Pj).sPlocal);
                end
            end
        case {'rev-rev'}
            Joints(Ji).mrows = 1; Joints(Ji).nbody = 2;
            Pi = Joints(Ji).iPindex; Pj = Joints(Ji).jPindex;
            Joints(Ji).iBindex = Points(Pi).Bindex;
            Joints(Ji).jBindex = Points(Pj).Bindex;
        case {'rev-tran'}
            Joints(Ji).mrows = 1; Joints(Ji).nbody = 2;
            Pi = Joints(Ji).iPindex; Pj = Joints(Ji).jPindex;
            Joints(Ji).iBindex = Points(Pi).Bindex;
            Joints(Ji).jBindex = Points(Pj).Bindex;
        case {'rel-rot', 'rel-tran'}
            Joints(Ji).mrows = 1; Joints(Ji).nbody = 1;
        case {'disc'}
            Joints(Ji).mrows = 2; Joints(Ji).nbody = 1;
        case {'rigid'}
            Joints(Ji).mrows = 3; Joints(Ji).nbody = 2;
            Bi = Joints(Ji).iBindex;    Bj = Joints(Ji).jBindex;
            if Bi == 0
                Joints(Ji).d0 = -Bodies(Bj).A'*Bodies(Bj).r;
                Joints(Ji).p0 = -Bodies(Bj).p;
            elseif Bj == 0
                Joints(Ji).d0 = Bodies(Bi).r;
                Joints(Ji).p0 = Bodies(Bi).p;
            else
                Joints(Ji).d0 = Bodies(Bj).A'*(Bodies(Bi).r - Bodies(Bj).r);
                Joints(Ji).p0 = Bodies(Bi).p - Bodies(Bj).p;
            end
            
       otherwise
            'Undefined joint type'
    end
end

% Functions
    nFc = length(Functs);
    for Ci = 1:nFc
        functData(Ci);
    end
% ------------------------------------------------------------
% Compute number of constraints and determine row/column pointers
    nConst = 0;
for Ji = 1:nJ
    Joints(Ji).rows = nConst + 1;
    Joints(Ji).rowe = nConst + Joints(Ji).mrows;
    nConst = Joints(Ji).rowe;
        Bi = Joints(Ji).iBindex;
        if Bi ~= 0
            Joints(Ji).colis = 3*(Bi - 1) + 1;
            Joints(Ji).colie = 3*Bi; 
        end
        Bj = Joints(Ji).jBindex;
        if Bj ~= 0
            Joints(Ji).coljs = 3*(Bj - 1) + 1;
            Joints(Ji).colje = 3*Bj; 
        end
end
