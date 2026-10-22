function CutJoint(check, theta, theta_d, t)
% Position constraints
    include_global

% Jacobian
    C = zeros(nConst, nJC);
    % cut 1: revolute C2 - C4
    C(1:2, 1) = C(1:2, 1) + [1; 0];
    C(1:2, 2) = C(1:2, 2) + [0; 1];
    C(1:2, 3) = C(1:2, 3) + s_rot(Points(7).rP - Bodies(1).r);
    C(1:2, 5) = C(1:2, 5) + s_rot(Points(7).rP - Points(2).rP);
    C(1:2, 4) = C(1:2, 4) + s_rot(Points(7).rP - Points(6).rP);
    C(1:2, 1) = C(1:2, 1) - [1; 0];
    C(1:2, 2) = C(1:2, 2) - [0; 1];
    C(1:2, 3) = C(1:2, 3) - s_rot(Points(8).rP - Bodies(1).r);
    C(1:2, 6) = C(1:2, 6) - s_rot(Points(8).rP - Points(4).rP);

if check == 1
    % Constraints
    Phi = zeros(nConst, 1);
    Phi(1:2) = Points(7).rP - Points(8).rP;

elseif check == 2
    % Right-hand-side of velocity constraints
    rhsV = zeros(nConst, 1);
    rhsV(1:2) = [0; 0];

elseif check == 3
    % Right-hand-side of acceleration constraints
    rhsA = zeros(nConst, 1);
    rhsA(1:2) = -(s_rot((Points(7).rP_d - Bodies(1).r_d)*theta_d(3) ...
            + (Points(7).rP_d - Points(2).rP_d)*theta_d(5) ...
            + (Points(7).rP_d - Points(6).rP_d)*theta_d(4)) ...
                - s_rot((Points(8).rP_d - Bodies(1).r_d)*theta_d(3) ...
            + (Points(8).rP_d - Points(4).rP_d)*theta_d(6)));
end

end
