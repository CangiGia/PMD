% Update_Position

% Compute A's
    for Bi=1:nB
        Bodies(Bi).A = Matrix_A(Bodies(Bi).p);
    end  
    
% Compute sP = A * sP_prime; rP = r + sP
    for Pi = 1:nP
        Bi = Points(Pi).Bindex;
        if Bi ~= 0
            Points(Pi).sP = Bodies(Bi).A * Points(Pi).sPlocal;
            Points(Pi).sP_r = s_rot(Points(Pi).sP);
            Points(Pi).rP = Bodies(Bi).r + Points(Pi).sP;
        end
    end
