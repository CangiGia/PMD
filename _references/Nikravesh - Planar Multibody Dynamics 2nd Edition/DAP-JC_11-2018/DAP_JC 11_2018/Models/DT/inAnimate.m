function inAnimate
    include_global

C1 = Point_struct;
C1.Bindex = 5;
C1.sPlocal = [0; 0];

Points_anim = [C1];

Bodies(1).color = 'r';
Bodies(2).color = 'g';
Bodies(3).color = 'b';
Bodies(4).color = 'm';
Bodies(5).color = 'c';

% Variables for defining the animation axes used by plot_system
xmin = -0.8; xmax = 2.2;
ymin = -0.4; ymax = 2;
