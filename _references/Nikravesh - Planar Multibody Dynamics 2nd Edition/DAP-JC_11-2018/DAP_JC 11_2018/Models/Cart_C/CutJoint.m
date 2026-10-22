function CutJoint(check, theta, theta_d, t)
% Position constarints
    include_global

% Jacobian
    R = 0.1; a = 0.6;
    C = [-R  0   R
          0  a   0];
      
if check == 1
    % Constraints
    Phi = Points(3).rP - Bodies(3).r;
    
elseif check == 2    
    % Right-hand-side of velocity constraints
    rhsV = [0; 0];
        
elseif check == 3         
    % Right-hand-side of acceleration constraints    
    rhsA = [0; 0];
end