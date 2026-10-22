function inBodies
    include_global

% Mountain bike 1 (Nikravesh, Section 15.9): body-coordinate model

B1 = Body_struct;
B1.m = 3.2;
B1.J = 0.1; % assumed: not given in the book
B1.r = [0.85; 0.55];
B1.p = 0;

B2 = Body_struct;
B2.m = 1.8;
B2.J = 0.03; % assumed: not given in the book
B2.r = [0.38; 0.33];
B2.p = 0;

B3 = Body_struct;
B3.m = 1.5;
B3.J = 0.16; % assumed (m R^2)
B3.r = [0.16; 0.33];
B3.p = 0;

B4 = Body_struct;
B4.m = 1.5;
B4.J = 0.16; % assumed (m R^2)
B4.r = [1.25399559; 0.33001102];
B4.p = 0;

Bodies = [B1; B2; B3; B4];
