function inBodies
    include_global

% Mountain bike 2 (Nikravesh, Section 15.9): body-coordinate model

B1 = Body_struct;
B1.m = 3.2;
B1.J = 0.1; % assumed: not given in the book
B1.r = [0.8; 0.54];
B1.p = 0;

B2 = Body_struct;
B2.m = 0.5;
B2.J = 0.008; % assumed (slender bar)
B2.r = [0.38; 0.29];
B2.p = 0;

B3 = Body_struct;
B3.m = 0.5;
B3.J = 0.0045; % assumed (slender bar)
B3.r = [0.221501; 0.44310985];
B3.p = 0;

B4 = Body_struct;
B4.m = 0.8;
B4.J = 0.0064; % assumed
B4.r = [0.45866785; 0.58432392];
B4.p = 0;

B5 = Body_struct;
B5.m = 1.2;
B5.J = 0.1; % assumed (m R^2)
B5.r = [0.16; 0.29];
B5.p = 0;

B6 = Body_struct;
B6.m = 1.2;
B6.J = 0.1; % assumed (m R^2)
B6.r = [1.24375; 0.29];
B6.p = 0;

Bodies = [B1; B2; B3; B4; B5; B6];
