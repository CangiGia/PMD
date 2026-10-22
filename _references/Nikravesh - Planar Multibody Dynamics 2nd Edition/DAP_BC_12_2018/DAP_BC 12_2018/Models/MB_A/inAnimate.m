function inAnimate
    include_global

C1 = Point_struct;
C1.Bindex = 1;
C1.sPlocal = [0.26; 0.14];

C2 = Point_struct;
C2.Bindex = 2;
C2.sPlocal = [0.14090909; 0.19453696];

Points_anim = [C1; C2];

Bodies(3).shape = 'circle';
Bodies(3).R = 0.33;
Bodies(4).shape = 'circle';
Bodies(4).R = 0.33;
Bodies(1).color = 'r';
Bodies(2).color = 'g';
Bodies(3).color = 'b';
Bodies(4).color = 'b';

% Variables for defining the animation axes used by plot_system
xmin = 0; xmax = 1.6;
ymin = 0; ymax = 1;
