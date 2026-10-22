function inAnimate
    include_global

    Points_anim = [];

Bodies(1).shape = 'circle';
Bodies(1).R = 0.05;

Bodies(2).shape = 'circle';
Bodies(2).R = 0.03;

Bodies(3).shape = 'circle';
Bodies(3).R = 0.03;

% Parameters for defining animation/plot axes 
    xmin = -0.5; xmax =  0.5;
    ymin = -0.5; ymax =  0.5;
