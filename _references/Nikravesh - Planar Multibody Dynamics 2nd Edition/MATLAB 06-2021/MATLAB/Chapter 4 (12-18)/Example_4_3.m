% Example 4.3
% Define coordinates
    r = [2.5; 1.2]; phi = 5.6723;
% Define constant coordinates
    s_A_local = [2.18; 0]; s_B_local = [-1.8; 1.3];
% Define mass and polar moment of inertia
    m = 2.5; J = 1.2;
% Define forces at A and B
    f_A = [2; -1]; f_B = [-3; 2];
% Compute matrix A
    A = A_matrix(phi);
% Compute global coordinates of A and B
    s_A = A*s_A_local; s_B = A*s_B_local;
% Compute the sum of forces
    f = f_A + f_B;
% Compute the sum of moments
    n = s_rot(s_A)'*f_A + s_rot(s_B)'*f_B;
% Compute accelerations
    r_dd = f/m
    phi_dd = n/J
    