function inBodies
    include_global

% Mountain bike 2 (Nikravesh, Section 15.9): joint-coordinate model

B1 = Body_struct;
B1.m = 3.2;
B1.J = 0.1; % assumed: not given in the book
B1.theta = [0.8; 0.54; 0];

B2 = Body_struct;
B2.m = 0.5;
B2.J = 0.008; % assumed (slender bar)
B2.theta = 0;

B3 = Body_struct;
B3.m = 0.5;
B3.J = 0.0045; % assumed (slender bar)
B3.theta = 0;

B4 = Body_struct;
B4.m = 0.8;
B4.J = 0.0064; % assumed
B4.theta = 0;

B5 = Body_struct;
B5.m = 1.2;
B5.J = 0.1; % assumed (m R^2)
B5.theta = 0;

B6 = Body_struct;
B6.m = 1.2;
B6.J = 0.1; % assumed (m R^2)
B6.theta = [1.24375; 0.29; 0];

Bodies = [B1; B2; B3; B4; B5; B6];
