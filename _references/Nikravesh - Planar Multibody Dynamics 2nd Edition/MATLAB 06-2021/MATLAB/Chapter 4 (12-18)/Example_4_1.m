% Example 4.1
% Define coordinates, masses and applied forces for the particles
    r1 = [2; 1];   r2 = [1; 4];   r3 = [-3; 2];
    m1 = 2;        m2 = 5;        m3 = 3;
    f1 = [-6; -2]; f2 = [0; 0];   f3 = [5; -1];    
% Compute the total mass
    m = m1 + m2 + m3;
% Compute the coordinates of mass center
    rC = (m1*r1 + m2*r2 + m3*r3)/m
% Compute the total force vector
    f = f1 + f2 + f3;
% Compute the acceleration of mass center
    rC_dd = f/m
% Compute s_i vectors
    s1 = r1 - rC; s2 = r2 - rC; s3 = r3 - rC;
% Compute the first moment equation
    sigma_m_s = m1*s1 + m2*s2 + m3*s3
    