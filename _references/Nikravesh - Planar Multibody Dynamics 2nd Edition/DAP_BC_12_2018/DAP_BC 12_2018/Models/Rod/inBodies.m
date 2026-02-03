function inBodies
    include_global

B1 = Body_struct;
B1.m = 1.0;
B1.J = 0.01;
B1.r = [0; 1];
B1.p = pi/4;
B1.r_d = [0; -6];

Bodies = [B1];