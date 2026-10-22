function inBodies
    include_global
    
B1 = Body_struct;
B1.m = 2;
B1.J = 0.5;
B1.p = 6.0819;

B2 = Body_struct;
B2.theta = 6.0819;
B2.m = 20;
B2.J = 2.5;

Bodies = [B1; B2];
