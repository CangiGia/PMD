function inAnimate
    include_global

    Points_anim = [];

Bodies(3).shape = 'rect';
Bodies(3).W = 0.04;
Bodies(3).H = 0.02;

% Bodies(2).shape = 'line';
% Bodies(2).H = 0.26;
% 
% Bodies(1).shape = 'line';
% Bodies(1).H = 0.12;

% Parameters for defining animation/plot axes 
    xmin = -0.3; xmax =  0.5;
    ymin = -0.4; ymax =  0.4;
