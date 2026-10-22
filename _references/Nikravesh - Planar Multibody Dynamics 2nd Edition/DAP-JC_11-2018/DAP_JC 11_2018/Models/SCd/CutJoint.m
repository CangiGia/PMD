function CutJoint(check, theta, theta_d, t)
% Position constraints
    include_global

% Jacobian
    dcut1 = Points(5).rP;
    dcut2 = Points(5).rP - Points(4).rP;
    u3 = Uvectors(2).u;
    C = [s_rot(dcut1)  s_rot(dcut2)   -u3
         1   0   0];

if check == 1
    % Constraints
    Phi = [(Points(5).rP - Points(6).rP)
           (theta(1) - pi/3 -2*pi*t)];
    
elseif check == 2    
    % Right-hand-side of velocity constraints
    rhsV = [0
            0
            2*pi];
        
elseif check == 3         
    % Right-hand-side of acceleration constraints    
    dcut1d = Points(5).rP_d;
    dcut2d = Points(5).rP_d - Points(4).rP_d;
    rhsA = -[s_rot(dcut1d*theta_d(1) + dcut2d*theta_d(2))
              0];
end

end