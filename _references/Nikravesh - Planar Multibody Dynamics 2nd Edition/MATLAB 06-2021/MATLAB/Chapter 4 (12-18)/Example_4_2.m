% Example 4.2
% Define coordinates
    r = [2.1; 1.6]; phi = pi/6;
% Define constant coordinates
    s_P_local = [-0.2; 0.3];
% Define mass 
    m = 0.2;
% Define force at P
    f_P = [1.2; 1.0]; T = -0.15;
% Compute matrix A
    A = A_matrix(phi);
% Compute global coordinates of A and B
    s_P = A*s_P_local; 
% Compute the sum of moments
    n = s_rot(s_A)'*f_P
% Compute weight
    weight = m*9.81*[0; -1]
% Construct array of force/moment
    h = [f_P + weight; n + T]
    