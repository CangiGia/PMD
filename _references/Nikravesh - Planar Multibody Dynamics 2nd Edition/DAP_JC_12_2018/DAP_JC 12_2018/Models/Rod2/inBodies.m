function inBodies
    include_global
    
B1 = Body_struct;
B1.m = 1.0;
B1.J = 0.1;
B1.theta = [0 1.0 0];

B2 = Body_struct;
B2.m = 1.0;
B2.J = 0.1;
B2.theta = -pi/4;

Bodies = [B1; B2];
