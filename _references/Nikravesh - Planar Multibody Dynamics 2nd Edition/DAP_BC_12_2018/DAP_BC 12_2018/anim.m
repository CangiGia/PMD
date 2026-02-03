% Animate response
    include_global
    
    nt = size(T,1);
    nP = nPtot;
    pause on
    
for Bi = 1:nB
    switch (Bodies(Bi).shape);
        case {'circle'}
            [cx cy cz] = cylinder(Bodies(Bi).R, 40);
            Bodies(Bi).circ = [cx(1,:); cy(1,:)];
        case {'rect'}
            w2 = Bodies(Bi).W/2; h2 = Bodies(Bi).H/2;
            Bodies(Bi).P4 = [w2 -w2 -w2  w2
                             h2  h2 -h2 -h2];
        case {'line'}
            if Bodies(Bi).W == 0
                h2 = Bodies(Bi).H/2;
                Bodies(Bi).P4 = [0  0
                                 h2 -h2];
            else
                w2 = Bodies(Bi).W/2;
                Bodies(Bi).P4 = [w2 -w2
                                 0   0];
            end
    end
end

for i = 1:nt
    u = uT(i,:)';
    u_to_Bodies; Update_Position;         
    plot_system,     drawnow;         
    if T(i) == 0 && nt > 1
        disp('  Press a key to continue!')
        pause 
    end
    pause off
end

