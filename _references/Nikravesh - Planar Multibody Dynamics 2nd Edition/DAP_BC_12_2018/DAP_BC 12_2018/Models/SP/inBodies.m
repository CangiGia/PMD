function inBodies
    include_global

B1 = Body_struct;
B1.m = 5.0;
B1.J = 4.0;
B1.r = [1.0; 0.2];

B2 = Body_struct;
B2.m = 2.0;
B2.J = 0.2;
B2.r = [1.25; -0.233];
B2.p = pi/6;

Bodies = [B1; B2];
