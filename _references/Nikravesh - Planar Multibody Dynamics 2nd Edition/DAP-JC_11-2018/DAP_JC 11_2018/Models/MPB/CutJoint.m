function CutJoint(check, theta, theta_d, t)
% Position constarints
    include_global

% Jacobian
    dcut1 = Points(6).rP - Points(1).rP;
    dcut2 = Points(6).rP - Points(3).rP;
    C  = Uvectors(1).u'*[dcut1  dcut2];if check == 1
    % Constraints
    Phi = s_rot(Uvectors(1).u)'*(Points(4).rP - Points(6).rP);
    
elseif check == 2    
    % Right-hand-side of velocity constraints
    rhsV = 0;
        
elseif check == 3         
    % Right-hand-side of acceleration constraints    
    dcut1 = Points(6).rP - Points(1).rP;
    dcut2 = Points(6).rP - Points(3).rP;
    
    dcut1d = Points(6).rP_d;
    dcut2d = Points(6).rP_d - Points(3).rP_d;
    
    rhsA = -(Uvectors(1).u_d'*[dcut1  dcut2] + ...
             Uvectors(1).u'*[dcut1d  dcut2d])*theta_d;
end