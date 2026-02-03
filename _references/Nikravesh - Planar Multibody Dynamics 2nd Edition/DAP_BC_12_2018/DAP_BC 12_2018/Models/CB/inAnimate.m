function inAnimate
    include_global

%-----------------------------------------------------------%
% Required data and the corresponding defauls values        %
%  Bi = Body_struct                                         %
%  Bi.shape = ' ';  circle, rect, line (no default shape)   %
%  Bi.W = 0;        width of rectangle or line              %
%  Bi.H = 0;        hight of rectangle or line              %
%  Bi.R = 0;        radius of circle/disc                   %
%                                                           %
% Output: None                                              %
%                                                           %
%  Pi = Body_struct                                         %
%  Pi.Bindex = 0;   body index                              %
%  Pi.sPlocal = [0;0];  body-fixed coordinates                  %
%                                                           %
% Output:                                                   %
%   Ppoints_anim = [     ]                                  %
%-----------------------------------------------------------%
    
Points_anim = [];

Bodies(1).shape = 'rect';
Bodies(1).W = 0.4;
Bodies(1).H = 0.2;


% Parameters for defining animation/plot axes 
    xmin = -0.5; xmax =  2.0;
    ymin = -1.0; ymax =  0.5;
