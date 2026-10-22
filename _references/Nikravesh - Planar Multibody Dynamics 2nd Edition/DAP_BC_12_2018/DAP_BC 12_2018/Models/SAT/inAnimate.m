function inAnimate
    include_global

Points_anim = [];

Bodies(1).shape = 'line';
Bodies(1).W = 0.6;
Bodies(1).color = 'r';
Bodies(2).color = 'g';
Bodies(3).color = 'b';
Bodies(4).color = 'm';
Bodies(5).color = 'c';

% Variables for defining the animation axes used by plot_system
xmin = -1.1; xmax = 0.2;
ymin = -0.5; ymax = 0.3;
