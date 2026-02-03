% Example 6.2 (two constrained bodies: continuation of Example 6.1)
    clc; clear all

%% Exaample 6.1
% Constant data
    m1 = 2.0; J1 = 1.5; m2 = 1.5; J2 = 0.5;
    k = 100; L0 = 2.0; dc = 20;
    g = 9.81;
% Coordinates and velocities (variables)
    rA1 = [-1.5; 3.0]; rB1  = [-1.5; 0]; 
    rC2 = [-1.5; 4.5]; rB1_d = [-2.0; 1.0];
    rG1 = [-1.5; 1.5]; rG2  = [-0.5; 4.5];
% Unit vectors along x and y axes 
    ux = [1; 0]; uy = [0; 1];
% Gravitational forces
    w1 = -g*m1*uy
    w2 = -g*m2*uy
% Spring force
    d1 = rC2 - rA1; L1 = norm(d1); u1 = d1/L1;
    fs = k*(L1 - L0)
    fsA =  fs*uy
    fsC = -fs*uy
% Damper force
    d2 = rB1; L2 = norm(d2); u2 = d2/L2;
    L2_d = u2'*rB1_d;
    fd = dc*L2_d
    fdB =  fd*ux
% Moment arms
    sA1 = rA1 - rG1;
    sB1 = rB1 - rG1;
    sC2 = rC2 - rG2;
% Array of forces for body 1 and body 2
    h1 = [(w1 + fsA + fdB)
          (s_rot(sA1)'*fsA + s_rot(sB1)'*fdB)]
    h2 = [(w2 + fsC)
          s_rot(sC2)'*fsC]
% Mass matrix
    M_diag = [m1 m1 J1 m2 m2 J2]';
    M = diag(M_diag)

%% Example 6.2
% Coordinates of D
    rD2 = [0.5; 4.5];
% massless link data
    Lc = 4.53;
% Unit vector along massless link
    uc = rD2/norm(rD2)
% Moment arm at D
    sD2 = rD2 - rG2;
    n2  = s_rot(sD2)'*uc
    
