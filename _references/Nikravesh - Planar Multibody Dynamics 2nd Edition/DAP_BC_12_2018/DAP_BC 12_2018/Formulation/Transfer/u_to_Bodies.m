% u_to_Bodies
% Unpack u into coordinate and velocity sub-arrays
    for Bi = 1:nB
        ir  = Bodies(Bi).irc;
        ird = Bodies(Bi).irv;
            Bodies(Bi).r   = u(ir:ir+1);
            Bodies(Bi).p   = u(ir+2);
            Bodies(Bi).r_d = u(ird:ird+1);
            Bodies(Bi).p_d = u(ird+2);
    end
