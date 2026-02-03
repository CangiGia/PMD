function [Phi, D] = Ex_12_4(theta)
% Constraints and Jacobian for a fourbar mechanism
    Phi = [(cos(theta(1)) + 3*cos(theta(2)) - 2.2*cos(theta(3)) - 2)
           (sin(theta(1)) + 3*sin(theta(2)) - 2.2*sin(theta(3)) - 0.5)
           (theta(1) - pi/2)];
    D   = [-sin(theta(1))  -3*sin(theta(2))  2.2*sin(theta(3))
            cos(theta(1))   3*cos(theta(2)) -2.2*cos(theta(3))
            1               0                0                ];
end
