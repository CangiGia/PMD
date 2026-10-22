function inBodies
    include_global

% Head and neck (Nikravesh, Section 15.14): joint-coordinate model

B1 = Body_struct;
B1.m = 50;
B1.J = 1.2;
B1.theta = 0;

B2 = Body_struct;
B2.m = 0.25;
B2.J = 0.0004;
B2.theta = 0;

B3 = Body_struct;
B3.m = 0.25;
B3.J = 0.0004;
B3.theta = 0;

B4 = Body_struct;
B4.m = 0.25;
B4.J = 0.0004;
B4.theta = 0;

B5 = Body_struct;
B5.m = 0.25;
B5.J = 0.0004;
B5.theta = 0;

B6 = Body_struct;
B6.m = 0.25;
B6.J = 0.0004;
B6.theta = 0;

B7 = Body_struct;
B7.m = 0.25;
B7.J = 0.0004;
B7.theta = 0;

B8 = Body_struct;
B8.m = 0.25;
B8.J = 0.0004;
B8.theta = 0;

B9 = Body_struct;
B9.m = 2;
B9.J = 0.03;
B9.theta = 0;

B10 = Body_struct;
B10.m = 1.4;
B10.J = 0.004;
B10.theta = [0.04610637; 0.251; 0];

Bodies = [B1; B2; B3; B4; B5; B6; B7; B8; B9; B10];
