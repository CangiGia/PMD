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

B3 = Body_struct;
B3.r = [0.4528; 0.6862];
B3.p = 5.0019;
B3.m = 0.5;
B3.J = 0.2;

Bodies = [B1; B2; B3];
