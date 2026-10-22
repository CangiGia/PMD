function CutJoint(check, theta, theta_d, t)
% Position constraints
    include_global

% Jacobian
    C = zeros(nConst, nJC);
    % cut 1: revolute Q4 - Q1
    C(1:2, 1) = C(1:2, 1) + [1; 0];
    C(1:2, 2) = C(1:2, 2) + [0; 1];
    C(1:2, 3) = C(1:2, 3) + s_rot(Points(9).rP - Bodies(1).r);
    C(1:2, 4) = C(1:2, 4) + s_rot(Points(9).rP - Points(2).rP);
    C(1:2, 5) = C(1:2, 5) + s_rot(Points(9).rP - Points(4).rP);
    C(1:2, 6) = C(1:2, 6) + s_rot(Points(9).rP - Points(8).rP);
    C(1:2, 1) = C(1:2, 1) - [1; 0];
    C(1:2, 2) = C(1:2, 2) - [0; 1];
    C(1:2, 3) = C(1:2, 3) - s_rot(Points(10).rP - Bodies(1).r);
    % cut 2: revolute-translational T1 - E6
    u2 = Uvectors(1).u; w2 = Uvectors(1).u_r;
    d2 = Points(13).rP - Points(14).rP;
    C(3, 1) = w2'*[1; 0];
    C(3, 2) = w2'*[0; 1];
    C(3, 3) = w2'*s_rot(Points(13).rP - Bodies(1).r) - u2'*d2;
    C(3, 8) = -w2'*[1; 0];
    C(3, 9) = -w2'*[0; 1];
    C(3, 10) = -w2'*s_rot(Points(14).rP - Bodies(6).r);

if check == 1
    % Constraints
    Phi = zeros(nConst, 1);
    Phi(1:2) = Points(9).rP - Points(10).rP;
    Phi(3) = w2'*d2 - 0;

elseif check == 2
    % Right-hand-side of velocity constraints
    rhsV = zeros(nConst, 1);
    rhsV(1:2) = [0; 0];
    rhsV(3) = 0;

elseif check == 3
    % Right-hand-side of acceleration constraints
    rhsA = zeros(nConst, 1);
    rhsA(1:2) = -(s_rot((Points(9).rP_d - Bodies(1).r_d)*theta_d(3) ...
            + (Points(9).rP_d - Points(2).rP_d)*theta_d(4) ...
            + (Points(9).rP_d - Points(4).rP_d)*theta_d(5) ...
            + (Points(9).rP_d - Points(8).rP_d)*theta_d(6)) ...
                - s_rot((Points(10).rP_d - Bodies(1).r_d)*theta_d(3)));
    dd2 = Points(13).rP_d - Points(14).rP_d;
    sa2 = s_rot((Points(13).rP_d - Bodies(1).r_d)*theta_d(3)) ...
          - s_rot((Points(14).rP_d - Bodies(6).r_d)*theta_d(10));
    rhsA(3) = 2*Bodies(1).p_d*(u2'*dd2) - w2'*sa2 + Bodies(1).p_d^2*(w2'*d2);
end

end
