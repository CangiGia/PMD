function CutJoint(check, theta, theta_d, t)
% Position constarints
    include_global

% Jacobian
    dcut1 = Points(4).rP - Points(1).rP;
    dcut2 = Points(4).rP - Points(3).rP;
    dcut3 = Points(4).rP - Points(6).rP;
    Cmat1 = [1 1 -1];
    Cmat2 = Uvectors(2).u'*[dcut1  dcut2   -dcut3];
    C = [Cmat1; Cmat2];
if check == 1
    % Constraints
    Phi = [s_rot(Uvectors(2).u)'*Uvectors(1).u
           s_rot(Uvectors(2).u)'*(Points(4).rP - Points(6).rP)];
    
elseif check == 2    
    % Right-hand-side of velocity constraints
    rhsV = [0
            0];
        
elseif check == 3         
    % Right-hand-side of acceleration constraints    
    dcut1d = Points(4).rP_d;
    dcut2d = Points(4).rP_d - Points(3).rP_d;
    dcut3d = Points(4).rP_d ;
    
    Cd2 = (Uvectors(2).u_d'*[dcut1  dcut2 -dcut3] + ...
           Uvectors(2).u'*[dcut1d  dcut2d -dcut3d])*theta_d;
    rhsA = -[0; Cd2];end

end