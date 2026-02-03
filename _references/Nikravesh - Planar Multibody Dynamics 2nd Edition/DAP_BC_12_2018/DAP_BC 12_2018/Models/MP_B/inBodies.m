function inBodies
    include_global

B1 = Body_struct;
B1.r = [0.5840; 0.3586];
B1.p = 6.0819;
B1.m = 20;
B1.J = 2.5;

B2 = Body_struct;
B2.m = 2;
B2.J = 0.5;
B2.r = [0.3450; 0.2900];
B2.p = 0;

Bodies = [B1; B2];
