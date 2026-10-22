function user_force
    include_global

% Shaking of the torso: horizontal force m*a*w^2*cos(w*t) (a = 0.1 m, f = 2 Hz)
% that moves the torso alone as x = x0 + a*(1 - cos(w*t)); t is read from Force_array
    t = evalin('caller', 't');
    a = 0.1; w = 2*pi*2;
    Bodies(1).f = Bodies(1).f + [Bodies(1).m*a*w^2*cos(w*t); 0];
