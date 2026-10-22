function inAnimate
    include_global

C1 = Point_struct;
C1.Bindex = 1;
C1.sPlocal = [-0.325; 0.02];

C2 = Point_struct;
C2.Bindex = 1;
C2.sPlocal = [0.351; 0.06];

Points_anim = [C1; C2];

Bodies(6).shape = 'circle';
Bodies(6).R = 0.31;
Bodies(7).shape = 'circle';
Bodies(7).R = 0.31;
Bodies(1).color = 'r';
Bodies(2).color = 'g';
Bodies(3).color = 'm';
Bodies(4).color = 'c';
Bodies(5).color = 'k';
Bodies(6).color = 'b';
Bodies(7).color = 'b';

% Variables for defining the animation axes used by plot_system
xmin = 0; xmax = 2.2;
ymin = 0; ymax = 1.2;
