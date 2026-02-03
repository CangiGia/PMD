% Example 11.2
clc; clear all

% Notation: m for minus or "before"
%           p for plus or "after"

% Define body-fixed vectors
    m1 = 0.2; m2 = 0.3;
    e = 0.6; u = [0; 1];
% Velocity vectors befor impact
    r1_d_m = [2; -2];
    r2_d_m = [1.5; 3];
% Equations
    DMD = [m1*eye(2) zeros(2) -u
           zeros(2) m2*eye(2)  u
           u'      -u'         0];
    rhs = [zeros(4,1);
         -(e + 1)*u'*(r1_d_m - r2_d_m)];
% Solve
    sol = DMD\rhs;
    r1_d_p = r1_d_m + sol(1:2)
    r2_d_p = r2_d_m + sol(3:4)
    impulse = sol(5)
% total energy before and after
    energy_m = (m1*r1_d_m'*r1_d_m + m2*r2_d_m'*r2_d_m)/2
    energy_p = (m1*r1_d_p'*r1_d_p + m2*r2_d_p'*r2_d_p)/2
    