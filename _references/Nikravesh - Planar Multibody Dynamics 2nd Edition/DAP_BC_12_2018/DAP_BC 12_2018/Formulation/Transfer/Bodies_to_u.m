% Bodies_to_u
% Pack coordinates and velocities ito u array
    for Bi = 1:nB
        ir  = Bodies(Bi).irc;
        ird = Bodies(Bi).irv;
            u(ir:ir+2)   = [Bodies(Bi).r
                            Bodies(Bi).p ];
            u(ird:ird+2) = [Bodies(Bi).r_d
                            Bodies(Bi).p_d ];
    end
