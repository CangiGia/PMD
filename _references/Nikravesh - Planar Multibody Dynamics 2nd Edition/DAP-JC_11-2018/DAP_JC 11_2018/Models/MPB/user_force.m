function user_force
    include_global

% Unilateral spring-damper in z-direction
    del = Bodies(2).r(2) - Forces(2).L0; % y - R
if del < 0
    fy = Forces(2).k*del + Forces(2).dc*Bodies(2).r_d(2);
    fsd = [0; -fy];
    Bodies(2).f = Bodies(2).f + fsd;
    Bodies(2).n = Bodies(2).n + s_rot(Points(5).sP)'*fsd;
end
