function inAnimate
    include_global

Points_anim = [];

Bodies(1).shape = 'rect';
Bodies(1).W = 0.4;
Bodies(1).H = 0.4;

Bodies(2).shape = 'line';
Bodies(2).H = 1.0;


% Parameters for defining animation/plot axes 
    xmin = -0.5; xmax =  2.0;
    ymin = -1.0; ymax =  0.5;
    