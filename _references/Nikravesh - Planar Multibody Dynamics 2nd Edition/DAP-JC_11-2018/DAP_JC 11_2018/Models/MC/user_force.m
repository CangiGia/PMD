function user_force
    include_global

% Unilateral spring-damper in the y-direction (tire of body 6)
    del = Bodies(6).r(2) - Forces(3).L0; % y-coordinate - radius
if del < 0
    fy = Forces(3).k*del + Forces(3).dc*Bodies(6).r_d(2);
    Bodies(6).f = Bodies(6).f + [0; -fy];
end

% Unilateral spring-damper in the y-direction (tire of body 7)
    del = Bodies(7).r(2) - Forces(3).L0; % y-coordinate - radius
if del < 0
    fy = Forces(3).k*del + Forces(3).dc*Bodies(7).r_d(2);
    Bodies(7).f = Bodies(7).f + [0; -fy];
end

