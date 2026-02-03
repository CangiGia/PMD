function Phi = s_c(x, theta1)
% Position constraints for a slider-crank mechanism
% x contains theta2 and theta3
global r_O L1 L3

% Unit vectors
    u1 = [cos(theta1); sin(theta1)];
    u3 = [cos(x(1)); sin(x(1))];
% Position constraints
    Phi = r_O + L1*u1 - L3*s_rot(u3) - x(2)*u3;
end

     