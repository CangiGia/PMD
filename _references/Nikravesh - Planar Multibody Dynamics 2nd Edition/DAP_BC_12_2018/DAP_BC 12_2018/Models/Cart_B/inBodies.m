function inBodies
    include_global

B1 = Body_struct;
B1.m = 20;
B1.J = 5;
B1.r = [0.5; 0.2];
B1.p = 0;

B2 = Body_struct;
B2.m = 2;
B2.J = 0.5;
B2.r = [0.2; 0.1];
B2.p = 0;

B3 = Body_struct;
B3.m = 2;
B3.J = 0.5;
B3.r = [0.8; 0.1];
B3.p = 0;

Bodies = [B1; B2; B3];
