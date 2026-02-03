% Example 4.4
% Define body-fixed vectors
    s_A1_local = [0.15; 0]; s_B2_local = [0; 0.1]; 
% Define coordinates (inertial) and velocities
    r1 = [ 0.1;  0.05]; phi1 = pi/4;
    r2 = [ 0.3; -0.05]; phi2 = 20*pi/180;
    r1_d = [ 0.1; 0.2]; phi1_d = -0.25;
    r2_d = [-0.2; 0.1]; phi2_d = 0.12;
% Define spring-damper-actiator constants
    L_0 = 0.15; k = 10; d_c = 5; f_a = -2;    
% Compute A matrices
    A1 = A_matrix(phi1); A2 = A_matrix(phi2);
% Compute global components of s_A1 and s_B2
    s_A1 = A1*s_A1_local; s_B2 = A2*s_B2_local;
% Compute global coordinates of A1 and B2
    r_A1 = r_Point(r1, s_A1);
    r_B2 = r_Point(r2, s_B2);
% Compute velocity of A1 and B2
    r_A1_d = r_Point_d(r1_d, s_A1, phi1_d);
    r_B2_d = r_Point_d(r2_d, s_B2, phi2_d);
% Compute vector d between A and B and its time derivative
    d   = r_A1 - r_B2;
    d_d = r_A1_d - r_B2_d;
% Compute force of the spring-damper 
    f_sda = pp_sda(d, d_d, k, L_0, d_c, f_a)
    f_A1 = -f_sda 
    f_B2 =  f_sda
% Compute each moment
    n1 = s_rot(s_A1)'*f_A1
    n2 = s_rot(s_B2)'*f_B2