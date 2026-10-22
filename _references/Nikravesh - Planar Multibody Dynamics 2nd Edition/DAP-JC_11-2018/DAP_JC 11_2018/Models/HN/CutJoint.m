function CutJoint(check, theta, theta_d, t)
% Position constraints
    include_global

% Jacobian
    C = zeros(nConst, nJC);

if check == 1
    % Constraints
    Phi = zeros(nConst, 1);


elseif check == 2
    % Right-hand-side of velocity constraints
    rhsV = zeros(nConst, 1);


elseif check == 3
    % Right-hand-side of acceleration constraints
    rhsA = zeros(nConst, 1);

end

end
