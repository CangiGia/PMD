function inBodies
    include_global
        
B1 = Body_struct;
B1.m = 2;
B1.J = 0.5;
B1.theta = 0;
B1.theta_d = 1.2;

B2 = Body_struct;
B2.m = 20;
B2.J = 5;
B2.theta = 0;

B3 = Body_struct;
B3.m = 2;
B3.J = 0.5;
B3.theta = 0;

Bodies = [B1; B2; B3];
