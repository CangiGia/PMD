% Example 3.1
% Define coordinates of the link
    r = [2.5; 1.2]; phi = 5.6723;
% Define constant coordinates (local)
    s_A_local = [2.18; 0]; s_B_local = [-1.8; 1.3];     
% Compute matrix A
    A = A_matrix(phi)
% Compute global components of s_A and s_B
    s_A = A*s_A_local
    s_B = A*s_B_local
% Compute global coordinates of A and B
    r_A = r_Point(r, s_A)
    r_B = r_Point(r, s_B)
% Compute s_BA
    s_BA_1 = A*(s_B_local - s_A_local)
    s_BA_2 = s_B - s_A
    s_BA_3 = r_B - r_A
    