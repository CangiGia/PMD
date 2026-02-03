% Example 11.4
    clc; clear all

% Define body-fixed vectors
    R = 0.1; m1 = 0.5; J1 = m1*R^2/2;
    L = 1; m2 = 0.6; J2 = m2*L^2/12;
    e = 0.8; u = [1; 0];
    sC1_local = [R; 0]; sC2_local = [0; L/4];

% Position and velocity befor impact
    c1_d_m = [2; 1; 0.5]; c2_d_m = [0; -1; 0];
    c_d_m = [c1_d_m; c2_d_m];
% Equations
    M1_array = [m1 m1 J1];
    M2_array = [m2 m2 J2];
    M = diag([M1_array M2_array]);
    A1 = eye(2); A2 = eye(2);
    sC1 = A1*sC1_local; sC2 = A2*sC2_local;

    D1 = [eye(2) s_rot(sC1)];
    D2 = [eye(2) s_rot(sC2)];
    Dij = [D1 -D2];
    DMD = [M      -Dij'*u
           u'*Dij       0];
    rhs = [zeros(6,1)
          -(e + 1)*u'*Dij*c_d_m];
% Solve
    sol = DMD\rhs;
    c1_d_p = c1_d_m + sol(1:3)
    c2_d_p = c2_d_m + sol(4:6)
    C1_d_p = D1*c1_d_p
    C2_d_p = D2*c2_d_p
    impulse = sol(7);
    