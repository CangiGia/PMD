function user_force
    include_global

% Unilateral spring-damper in the y-direction (tire of body 3)
    del = Bodies(3).r(2) - Forces(3).L0; % y-coordinate - radius
if del < 0
    fy = Forces(3).k*del + Forces(3).dc*Bodies(3).r_d(2);
    Bodies(3).f = Bodies(3).f + [0; -fy];
end

% Unilateral spring-damper in the y-direction (tire of body 4)
    del = Bodies(4).r(2) - Forces(3).L0; % y-coordinate - radius
if del < 0
    fy = Forces(3).k*del + Forces(3).dc*Bodies(4).r_d(2);
    Bodies(4).f = Bodies(4).f + [0; -fy];
end

