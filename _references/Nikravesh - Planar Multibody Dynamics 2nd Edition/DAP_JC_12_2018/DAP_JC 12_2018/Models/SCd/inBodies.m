function inBodies
    include_global
    
B1 = Body_struct;
B1.m = 1.0;
B1.J = 0.1;
B1.theta = pi/3;
B1.theta_d = 2*pi;

B2 = Body_struct;
B2.m = 2.0;
B2.J = 0.2;
B2.theta = -0.6;

B3 = Body_struct;
B3.m = 4.0;
B3.J = 0.4;
B3.theta = 0.4;

Bodies = [B1; B2; B3];
