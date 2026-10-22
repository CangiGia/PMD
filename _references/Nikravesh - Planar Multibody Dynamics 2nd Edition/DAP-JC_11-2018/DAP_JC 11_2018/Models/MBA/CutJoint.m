function CutJoint(check, theta, theta_d, t)
% Position constraints
    include_global

% Jacobian
    C = zeros(nConst, nJC);
    % cut 1: revolute-translational T1 - D4
    u1 = Uvectors(1).u; w1 = Uvectors(1).u_r;
    d1 = Points(7).rP - Points(8).rP;
    C(1, 1) = w1'*[1; 0];
    C(1, 2) = w1'*[0; 1];
    C(1, 3) = w1'*s_rot(Points(7).rP - Bodies(1).r) - u1'*d1;
    C(1, 6) = -w1'*[1; 0];
    C(1, 7) = -w1'*[0; 1];
    C(1, 8) = -w1'*s_rot(Points(8).rP - Bodies(4).r);

if check == 1
    % Constraints
    Phi = zeros(nConst, 1);
    Phi(1) = w1'*d1 - 0;

elseif check == 2
    % Right-hand-side of velocity constraints
    rhsV = zeros(nConst, 1);
    rhsV(1) = 0;

elseif check == 3
    % Right-hand-side of acceleration constraints
    rhsA = zeros(nConst, 1);
    dd1 = Points(7).rP_d - Points(8).rP_d;
    sa1 = s_rot((Points(7).rP_d - Bodies(1).r_d)*theta_d(3)) ...
          - s_rot((Points(8).rP_d - Bodies(4).r_d)*theta_d(8));
    rhsA(1) = 2*Bodies(1).p_d*(u1'*dd1) - w1'*sa1 + Bodies(1).p_d^2*(w1'*d1);
end

end
