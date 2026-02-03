function inAnimate
    include_global
    
C2 = Point_struct;
C2.Bindex = 2;
C2.sPlocal = [0.03; 0];

Points_anim = [C2];

Bodies(1).color = 'r';
Bodies(3).color = 'b';

% Variables for defining the 3D animation axes used by plot_sys

xmin =  0; xmax =  0.8;
ymin =  0; ymax =  0.8;
