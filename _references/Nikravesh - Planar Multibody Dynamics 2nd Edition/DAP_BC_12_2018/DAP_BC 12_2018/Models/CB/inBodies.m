function inBodies
    include_global

B1 = Body_struct;
B1.m = 1.0;
B1.J = 1.0;
B1.r = [1.0; 0.2];
B1.r_d = [ 0.0; 0.0];

Bodies = [B1];
