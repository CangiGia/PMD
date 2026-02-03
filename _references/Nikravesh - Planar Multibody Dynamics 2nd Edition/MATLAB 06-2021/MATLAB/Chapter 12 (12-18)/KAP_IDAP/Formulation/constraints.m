function [Phi, D] = constraints(theta)
    global L RQ theta1

% unit vectors and rotated unit vectors
    [u, ur] = uVectors(theta); 
% constraints
    Phi = [L(1)*u{1} + L(2)*u{2} - L(3)*u{3} - RQ
           theta(1) - theta1];
% Jacobian matrix
    D  = [L(1)*ur{1} L(2)*ur{2} -L(3)*ur{3}
          1            0             0        ]; 
