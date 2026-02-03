function ud = Ex_13_9_1(t, u) 
% Two balls: Equations of motion
    global M h_a L

% Jacobian and gamma
    p1 = u(3); u1 = [cos(p1); sin(p1)]; u1r = s_rot(u1);
    p2 = u(6); u2 = [cos(p2); sin(p2)]; u2r = s_rot(u2);
    D = [eye(2) -L*u1 zeros(2,3)
         zeros(2,3) -eye(2) L*u2];
    gamma1 =  L*u1r*u(9)^2;
    gamma2 = -L*u2r*u(12)^2;
    gamma  = [gamma1; gamma2];
% EQM
    DMD = [M -D'
           D  zeros(4)];
    rhs = [h_a; gamma];
% Compute  accelerations
    sol = DMD\rhs;
% Construct ud array
    ud = [u(7:12); sol(1:6)];
end
