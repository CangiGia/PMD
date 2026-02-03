% Bodies_to_u_d
% pack velocities and accelerations into u_dot
    u_d = zeros (nB6,1);
    for Bi = 1:nB
        ir  = Bodies(Bi).irc;
        ird = Bodies(Bi).irv;
            u_d(ir:ir+2)   = [Bodies(Bi).r_d
                              Bodies(Bi).p_d ];
            u_d(ird:ird+2) = [Bodies(Bi).r_dd
                              Bodies(Bi).p_dd];
    end
