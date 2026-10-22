function CutJoint(check, theta, theta_d, t)
% Position constarints
    include_global

% Jacobian
    dcut1 = Points(4).rP - Points(1).rP;
    dcut2 = Points(4).rP - Points(3).rP;
    dcut3 = Points(5).rP - Points(6).rP;
    C = [s_rot(dcut1)  s_rot(dcut2)   -s_rot(dcut3)];

if check == 1
    % Constraints
    Phi = Points(4).rP - Points(5).rP;
    
elseif check == 2    
    % Right-hand-side of velocity constraints
    rhsV = [0
            0];
        
elseif check == 3         
    % Right-hand-side of acceleration constraints    
    dcut1d = Points(4).rP_d;
    dcut2d = Points(4).rP_d - Points(3).rP_d;
    dcut3d = Points(5).rP_d ;
    rhsA = -s_rot(dcut1d*theta_d(1) + dcut2d*theta_d(2) ...
                                    - dcut3d*theta_d(3));
end

end