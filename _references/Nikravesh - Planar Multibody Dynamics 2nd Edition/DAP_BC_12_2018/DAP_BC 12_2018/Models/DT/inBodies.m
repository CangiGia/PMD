function inBodies
    include_global

% Dump truck (Nikravesh, Section 15.11): body-coordinate model

B1 = Body_struct;
B1.m = 3;
B1.J = 0.04;
B1.r = [-0.19696155; 0.03472964];
B1.p = 0;

B2 = Body_struct;
B2.m = 15;
B2.J = 5;
B2.r = [0.60141715; 0.16588428];
B2.p = 0;

B3 = Body_struct;
B3.m = 15;
B3.J = 5;
B3.r = [0.60002747; 0.19258825];
B3.p = 0;

B4 = Body_struct;
B4.m = 3;
B4.J = 0.04;
B4.r = [-0.20147324; 0.20985313];
B4.p = 0;

B5 = Body_struct;
B5.m = 300;
B5.J = 125;
B5.r = [0.78993311; 0.64835923];
B5.p = 0.01736307;

Bodies = [B1; B2; B3; B4; B5];
