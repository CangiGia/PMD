function inBodies
    include_global
    
B1 = Body_struct;
B1.m = 4.0;
B1.J = 0.8;
B1.theta = pi/6;
B1.theta_d = 0.3;

B2 = Body_struct;
B2.m = 1.5;
B2.J = 0.15;
B2.theta = -75*pi/180;
B2.theta_d = 0.1;

B3 = Body_struct;
B3.m = 1.5;
B3.J = 0.15;
B3.theta = pi/2;
B3.theta_d = -0.2;

Bodies = [B1; B2; B3];
