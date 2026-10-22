function inAnimate
    include_global

C1 = Point_struct;
C1.Bindex = 10;
C1.sPlocal = [0; 0];

Points_anim = [C1];

Bodies(1).color = 'r';
Bodies(9).color = 'b';
Bodies(10).color = 'm';

% Variables for defining the animation axes used by plot_system
xmin = -0.3; xmax = 0.3;
ymin = 0; ymax = 0.3;
