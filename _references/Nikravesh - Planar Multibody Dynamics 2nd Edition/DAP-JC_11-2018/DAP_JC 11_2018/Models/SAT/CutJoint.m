function CutJoint(check, theta, theta_d, t)
% Position constraints
    include_global

% Jacobian
    C = zeros(nConst, nJC);
    % cut 1: revolute C3 - C0
    C(1:2, 1) = C(1:2, 1) + s_rot(Points(10).rP - Points(4).rP);
    C(1:2, 2) = C(1:2, 2) + s_rot(Points(10).rP - Points(6).rP);
    C(1:2, 3) = C(1:2, 3) + s_rot(Points(10).rP - Points(8).rP);
    % cut 2: revolute D5 - D0
    C(3:4, 1) = C(3:4, 1) + s_rot(Points(14).rP - Points(4).rP);
    C(3:4, 2) = C(3:4, 2) + s_rot(Points(14).rP - Points(6).rP);
    C(3:4, 3) = C(3:4, 3) + s_rot(Points(14).rP - Points(8).rP);
    C(3:4, 4) = C(3:4, 4) + s_rot(Points(14).rP - Points(11).rP);
    C(3:4, 5) = C(3:4, 5) + s_rot(Points(14).rP - Points(13).rP);

if check == 1
    % Constraints
    Phi = zeros(nConst, 1);
    Phi(1:2) = Points(10).rP - Points(2).rP;
    Phi(3:4) = Points(14).rP - Points(3).rP;

elseif check == 2
    % Right-hand-side of velocity constraints
    rhsV = zeros(nConst, 1);
    rhsV(1:2) = [0; 0];
    rhsV(3:4) = [0; 0];

elseif check == 3
    % Right-hand-side of acceleration constraints
    rhsA = zeros(nConst, 1);
    rhsA(1:2) = -(s_rot((Points(10).rP_d - Points(4).rP_d)*theta_d(1) ...
            + (Points(10).rP_d - Points(6).rP_d)*theta_d(2) ...
            + (Points(10).rP_d - Points(8).rP_d)*theta_d(3)) ...
                - [0; 0]);
    rhsA(3:4) = -(s_rot((Points(14).rP_d - Points(4).rP_d)*theta_d(1) ...
            + (Points(14).rP_d - Points(6).rP_d)*theta_d(2) ...
            + (Points(14).rP_d - Points(8).rP_d)*theta_d(3) ...
            + (Points(14).rP_d - Points(11).rP_d)*theta_d(4) ...
            + (Points(14).rP_d - Points(13).rP_d)*theta_d(5)) ...
                - [0; 0]);
end

end
