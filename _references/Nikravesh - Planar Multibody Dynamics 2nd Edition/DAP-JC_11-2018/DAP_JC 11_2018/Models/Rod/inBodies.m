function inBodies
    include_global
    
B1 = Body_struct;
B1.m = 1.0;
B1.J = 0.01;
B1.theta = [0; 1; pi/4];
B1.theta_d = [0; -6; 0];

Bodies = [B1];