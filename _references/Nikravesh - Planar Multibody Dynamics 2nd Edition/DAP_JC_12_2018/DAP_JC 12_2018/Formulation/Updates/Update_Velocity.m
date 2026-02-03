% Update_Velocity

% Update body velocities
    for Bi = 1:nB
        irs = Bodies(Bi).irc;
        Bodies(Bi).r_d = c_d(irs:irs+1);
        Bodies(Bi).p_d = c_d(irs+2);
    end
% Compute sP_dot and rP_dot vectors
    for Pi = 1:nP
        Bi = Points(Pi).Bindex;
        if Bi ~= 0
            Points(Pi).sP_d = Points(Pi).sP_r * Bodies(Bi).p_d;
            Points(Pi).rP_d = Bodies(Bi).r_d + Points(Pi).sP_d;
        end
    end
        
% Compute u_dot vectors
    for Vi = 1:nU
        Bi = Uvectors(Vi).Bindex;
        if Bi ~= 0
            Uvectors(Vi).u_d = Uvectors(Vi).u_r * Bodies(Bi).p_d;
        end
    end
