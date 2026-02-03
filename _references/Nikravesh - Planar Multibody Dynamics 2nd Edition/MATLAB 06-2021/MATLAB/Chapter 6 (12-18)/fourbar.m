function Phi = fourbar(x, theta1)
% Position constraints for a four-bar mechanism
% Array x contains theta2 and theta3
global L0 L1 L2 L3

% Unit vectors
    u1 = [cos(theta1); sin(theta1)];
    u2 = [cos(x(1)); sin(x(1))];
    u3 = [cos(x(2)); sin(x(2))];
% Position constraints
    Phi = L1*u1 + L2*u2 - L3*u3 - [L0; 0];
end
