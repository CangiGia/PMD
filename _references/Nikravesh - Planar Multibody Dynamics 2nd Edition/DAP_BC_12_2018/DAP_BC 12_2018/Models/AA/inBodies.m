function inBodies
    include_global

B1 = Body_struct;
B1.m = 2;
B1.J = 0.5;
B1.r = [0.4398; 0.2512]; % correct
B1.p = -0.0367; % correct
% B1.r = [0.51; 0.28]; % incorrect
% B1.p = 340*pi/180; % incorrect

B2 = Body_struct;
B2.r = [0.6817; 0.3498]; % correct
B2.p = 0.0783; % correct
% B2.r = [0.75; 0.35]; % incorrect
% B2.p = 0; % incorrect
B2.m = 30;
B2.J = 2.5;

B3 = Body_struct;
B3.r = [0.4463; 0.4308]; % correct
B3.p = 6.5222; % correct
% B3.r = [0.49; 0.41]; % incorrect
% B3.p = 350*pi/180; % incorrect
B3.m = 1;
B3.J = 0.5;

Bodies = [B1; B2; B3];
