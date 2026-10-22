function inBodies
    include_global

% Dump truck (Nikravesh, Section 15.11): joint-coordinate model

B1 = Body_struct;
B1.m = 3;
B1.J = 0.04;
B1.theta = 0;

B2 = Body_struct;
B2.m = 15;
B2.J = 5;
B2.theta = 0;

B3 = Body_struct;
B3.m = 15;
B3.J = 5;
B3.theta = 0;

B4 = Body_struct;
B4.m = 3;
B4.J = 0.04;
B4.theta = 0;

B5 = Body_struct;
B5.m = 300;
B5.J = 125;
B5.theta = 0.01736307;

Bodies = [B1; B2; B3; B4; B5];
