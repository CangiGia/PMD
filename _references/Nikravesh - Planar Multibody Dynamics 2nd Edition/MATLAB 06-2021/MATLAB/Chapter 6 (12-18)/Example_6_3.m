% Example 6.3 (sliding pendulum)
    clc; clear all

% Constant data
    m1 = 5.0; J1 = 4.0; m2 = 2.0; J2 = 0.2;
    k = 20; L0 = 0.6;
    g = 9.81;
% Coordinates and velocities (variables)
    rO = [ 0; 0.2];  
    rG1 = [ 1.0; 0.2]; rG2  = [ 1.25;-0.233];
% Unit vectors along x and y axes 
    ux = [1; 0]; uy = [0; 1];
% Gravitational forces
    w1 = -g*m1*uy
    w2 = -g*m2*uy
% Spring force
    d1 = rG1 - rO; L1 = norm(d1); u1 = d1/L1;
    fs = k*(L1 - L0)
    fsG1 = -fs*u1
% Array of forces for body 1 and body 2
    h1 = [(w1 + fsG1)
           0]
    h2 = [ w2
           0]
% Mass matrix
    M_diag = [m1 m1 J1 m2 m2 J2]';
    M = diag(M_diag)

% Moment arm at B
    sB2 = rG1 - rG2;
    n2  = s_rot(sB2)'*uc
    