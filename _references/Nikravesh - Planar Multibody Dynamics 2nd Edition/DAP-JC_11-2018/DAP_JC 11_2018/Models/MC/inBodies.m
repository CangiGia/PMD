function inBodies
    include_global

% Motorcycle (Nikravesh, Section 15.10): joint-coordinate model

B1 = Body_struct;
B1.m = 230;
B1.J = 10; % assumed: not given in the book (mass center at (0.30, 0.15) in the frame datum A)
B1.theta = [0.98166089; 0.7014396; -0.31398685];

B2 = Body_struct;
B2.m = 10;
B2.J = 0.31; % assumed (slender bar)
B2.theta = 0;

B3 = Body_struct;
B3.m = 2;
B3.J = 0.0065; % assumed (slender bar)
B3.theta = 0.31398685;

B4 = Body_struct;
B4.m = 2;
B4.J = 0.004; % assumed (slender bar)
B4.theta = 0.31398685;

B5 = Body_struct;
B5.m = 8;
B5.J = 0.2; % assumed (slender bar)
B5.theta = 0.31398685;

B6 = Body_struct;
B6.m = 10;
B6.J = 0.96; % assumed (m R^2)
B6.theta = 0;

B7 = Body_struct;
B7.m = 10;
B7.J = 0.96; % assumed (m R^2)
B7.theta = 0;

Bodies = [B1; B2; B3; B4; B5; B6; B7];
