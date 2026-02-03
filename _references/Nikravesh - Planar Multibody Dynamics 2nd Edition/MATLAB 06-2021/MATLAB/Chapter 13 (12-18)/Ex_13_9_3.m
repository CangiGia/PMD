function vp = Ex_13_9_3(u, M, R, L, e)
% Perfome momentum balance: constrained multibody system

    vm = u(7:12); % velocities before impact
    p1 = u(3); u1 = [cos(p1); sin(p1)]; u1r = s_rot(u1);
    p2 = u(6); u2 = [cos(p2); sin(p2)]; u2r = s_rot(u2);
    D = [eye(2) -L*u1 zeros(2,3)
         zeros(2,3) -eye(2) L*u2];
    d12 = u(1:2) - u(4:5);
    u12 = d12/norm(d12);
    Du12 = u12'*[eye(2) -R*u1r -eye(2) -R*u2r];
    DD = [D; Du12];
    DMD = [M   -DD'
           DD   zeros(5,5)];
    rhs = [zeros(10,1); -(e + 1)*Du12*vm];
    sol = DMD\rhs;
    vp = vm + sol(1:6); % velocities after impact
end
