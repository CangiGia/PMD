% Example 11.3
    clc; clear all

% Data 
    L = 2; m = 1.0; J = 0.01;
    sC_local = [L/2; 0];
    e = 0.95; u = [0; 1];

% Position and velocity befor impact
    phi = 315*pi/180; 
    r = [0; -(L/2)*sin(phi)];
    c_d_m = [0; -6.5; 0];
% Equations
    M_array = [m m J];
    M = diag(M_array);
    A = A_matrix(phi);
    sC = A*sC_local;
    Dij = [eye(2) s_rot(sC)];
    DMD = [M      -Dij'*u
           u'*Dij    0];
    rhs = [zeros(3,1);
          -(e + 1)*u'*Dij*c_d_m];
% Solve
    sol = DMD\rhs;
    c_d_p = c_d_m + sol(1:3)
    impulse = sol(4);
    rC_d_p = Dij*c_d_p 
% total energy before and after
%     energy_m = (c_d_m'*M*c_d_m)/2
%     energy_p = (c_d_p'*M*c_d_p)/2
    