function inBodies
    include_global
    
B1 = Body_struct;
B1.m = 2;
B1.J = 0.5;
B1.theta = -0.0367; % correct
% B1.theta = [0.51; 0.28; 340*pi/180]; % incorrect

B2 = Body_struct;
B2.theta = [0.0783 + 0.0367]; % correct
% B2.theta = [0.0367]; % incorrect
B2.m = 30;
B2.J = 2.5;

B3 = Body_struct;
B3.theta = 6.5222; % correct
% B3.theta = 350*pi/180; % incorrect
B3.m = 1;
B3.J = 0.5;

Bodies = [B1; B2; B3];
