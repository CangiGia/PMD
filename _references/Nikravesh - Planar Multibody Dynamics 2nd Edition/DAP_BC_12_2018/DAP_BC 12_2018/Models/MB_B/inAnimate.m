function inAnimate
    include_global

C1 = Point_struct;
C1.Bindex = 1;
C1.sPlocal = [0.32; 0.08];

C2 = Point_struct;
C2.Bindex = 1;
C2.sPlocal = [-0.18; -0.17];

Points_anim = [C1; C2];

Bodies(5).shape = 'circle';
Bodies(5).R = 0.29;
Bodies(6).shape = 'circle';
Bodies(6).R = 0.29;
Bodies(1).color = 'r';
Bodies(2).color = 'g';
Bodies(3).color = 'm';
Bodies(4).color = 'c';
Bodies(5).color = 'b';
Bodies(6).color = 'b';

% Variables for defining the animation axes used by plot_system
xmin = 0; xmax = 1.6;
ymin = 0; ymax = 1;
