% Example 3.2
% Define coordinates
    r = [2.5; 1.2]; phi = 5.6723;
% Define constant coordinates (local)
     s_A_local = [2.18; 0]; s_B_local = [-1.8; 1.3];
% Define velocities and accelerations
    r_d = [1; -2]; phi_d = 1;
    r_dd = [-0.4; 0.4]; phi_dd = 4.65;    
% Compute matrix A
    A = A_matrix(phi);
% Compute global components of s_A and s_B
    s_A = A*s_A_local; s_B = A*s_B_local;
% Compute global coordinates of A and B
    r_A = r_Point(r, s_A);
    r_B = r_Point(r, s_B);
% Compute veloccity of A and B
    r_A_d = r_Point_d(r_d, s_A, phi_d)
    r_B_d = r_Point_d(r_d, s_B, phi_d)
% Compute acceleration of A and B
    r_A_dd = r_Point_dd(r_dd, s_A, phi_d, phi_dd)
    r_B_dd = r_Point_dd(r_dd, s_B, phi_d, phi_dd)