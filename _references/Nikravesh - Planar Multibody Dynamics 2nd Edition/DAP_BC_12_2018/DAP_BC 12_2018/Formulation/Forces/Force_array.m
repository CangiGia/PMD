function g = Force_array(t)
    include_global    

% initialize body force vectors
for Bi=1:nB
    Bodies(Bi).f = [0; 0]; Bodies(Bi).n = 0;
end
for Fi = 1:nF
    switch (Forces(Fi).type);
        case {'weight'}
            for Bi=1:nB
                Bodies(Bi).f = Bodies(Bi).f + Bodies(Bi).wgt;
            end
        case {'ptp'}
            SDA_ptp
        case {'rot-sda'}
            SDA_rot
        case {'flocal'}
            Bi = Forces(Fi).iBindex;
            Bodies(Bi).f = Bodies(Bi).f + Bodies(Bi).A*Forces(Fi).flocal; 
        case {'f'}
            Bi = Forces(Fi).iBindex;
            Bodies(Bi).f = Bodies(Bi).f + Forces(Fi).f; 
        case {'T'}
            Bi = Forces(Fi).iBindex;
            Bodies(Bi).n = Bodies(Bi).n + Forces(Fi).T; 
        case {'user'}
            user_force
    end
end

    g = zeros(nB3,1);
for Bi = 1:nB
    ks = Bodies(Bi).irc; ke = ks + 2;
    g(ks:ke) = [Bodies(Bi).f; Bodies(Bi).n];   
end
