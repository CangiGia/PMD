function user_force
    include_global

% Unilateral spring-damper in the y-direction
    del = Bodies(2).r(2) - Forces(2).L0; % y-coordinate - Radius
if del < 0
    fy = Forces(2).k*del + Forces(2).dc*Bodies(2).r_d(2);
    fsd = [0; -fy];
    Bodies(2).f = Bodies(2).f + fsd;
end
