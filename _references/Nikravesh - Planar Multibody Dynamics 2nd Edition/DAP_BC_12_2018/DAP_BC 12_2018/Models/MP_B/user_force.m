function user_force
    include_global

% Unilateral spring-damper in z-direction
    del = Bodies(1).r(2) - Forces(2).L0; % y - R
if del < 0
    fy = Forces(2).k*del + Forces(2).dc*Bodies(1).r_d(2);
    fsd = [0; -fy];
    Bodies(1).f = Bodies(1).f + fsd;
    Bodies(1).n = Bodies(1).n + s_rot(Points(3).sP)'*fsd;
end
