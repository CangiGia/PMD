function CutJoint(check, theta, theta_d, t)
% Position constraints
    include_global

% Jacobian
    C = zeros(nConst, nJC);
    % cut 1: revolute D2 - D3
    C(1:2, 1) = C(1:2, 1) + s_rot(Points(6).rP - Points(3).rP);
    C(1:2, 2) = C(1:2, 2) + s_rot(Points(6).rP - Points(5).rP);
    C(1:2, 3) = C(1:2, 3) - s_rot(Points(8).rP - Points(10).rP);
    % cut 2: revolute K2 - K5
    C(3:4, 1) = C(3:4, 1) + s_rot(Points(7).rP - Points(3).rP);
    C(3:4, 2) = C(3:4, 2) + s_rot(Points(7).rP - Points(5).rP);
    C(3:4, 3) = C(3:4, 3) - s_rot(Points(14).rP - Points(10).rP);
    C(3:4, 4) = C(3:4, 4) - s_rot(Points(14).rP - Points(11).rP);
    C(3:4, 5) = C(3:4, 5) - s_rot(Points(14).rP - Points(13).rP);
    % cut 3: relative-translational driver D2 - A0
    d3 = Points(6).rP - Points(1).rP;
    [f3, f3_d, f3_dd] = functs(1, t);
    C(5, 1) = C(5, 1) + d3'*s_rot(Points(6).rP - Points(3).rP);
    C(5, 2) = C(5, 2) + d3'*s_rot(Points(6).rP - Points(5).rP);

if check == 1
    % Constraints
    Phi = zeros(nConst, 1);
    Phi(1:2) = Points(6).rP - Points(8).rP;
    Phi(3:4) = Points(7).rP - Points(14).rP;
    Phi(5) = (d3'*d3 - f3^2)/2;

elseif check == 2
    % Right-hand-side of velocity constraints
    rhsV = zeros(nConst, 1);
    rhsV(1:2) = [0; 0];
    rhsV(3:4) = [0; 0];
    rhsV(5) = f3*f3_d;

elseif check == 3
    % Right-hand-side of acceleration constraints
    rhsA = zeros(nConst, 1);
    rhsA(1:2) = -(s_rot((Points(6).rP_d - Points(3).rP_d)*theta_d(1) ...
            + (Points(6).rP_d - Points(5).rP_d)*theta_d(2)) ...
                - s_rot((Points(8).rP_d - Points(10).rP_d)*theta_d(3)));
    rhsA(3:4) = -(s_rot((Points(7).rP_d - Points(3).rP_d)*theta_d(1) ...
            + (Points(7).rP_d - Points(5).rP_d)*theta_d(2)) ...
                - s_rot((Points(14).rP_d - Points(10).rP_d)*theta_d(3) ...
            + (Points(14).rP_d - Points(11).rP_d)*theta_d(4) ...
            + (Points(14).rP_d - Points(13).rP_d)*theta_d(5)));
    dd3 = Points(6).rP_d - Points(1).rP_d;
    sa3 = s_rot((Points(6).rP_d - Points(3).rP_d)*theta_d(1) ...
            + (Points(6).rP_d - Points(5).rP_d)*theta_d(2)) ...
          - [0; 0];
    rhsA(5) = f3_d^2 + f3*f3_dd - dd3'*dd3 - d3'*sa3;
end

end
