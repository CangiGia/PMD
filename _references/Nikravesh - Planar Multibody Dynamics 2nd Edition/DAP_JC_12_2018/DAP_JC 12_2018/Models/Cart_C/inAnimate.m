function inAnimate
    include_global
    
P5 = Point_struct;
P5.Bindex = 1;
P5.sPlocal = [0.1; 0];

P6 = Point_struct;
P6.Bindex = 3;
P6.sPlocal = [0.1; 0];

Points_anim = [P5; P6];

Bodies(2).shape = 'rect';
Bodies(2).W = 0.9;
Bodies(2).H = 0.2;
Bodies(2).color = 'r';

Bodies(1).shape = 'circle';
Bodies(1).R = 0.1;

Bodies(3).shape = 'circle';
Bodies(3).R = 0.1;

% Variables for defining the 3D animation axes used by plot_sys

xmin = 0; xmax =  1.5;
ymin = -0.5; ymax =  1.0;
