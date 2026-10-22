function inBodies
    include_global

% Motorcycle (Nikravesh, Section 15.10): body-coordinate model

B1 = Body_struct;
B1.m = 230;
B1.J = 10; % assumed: not given in the book (mass center at (0.30, 0.15) in the frame datum A)
B1.r = [0.98166089; 0.7014396];
B1.p = -0.31398685;

B2 = Body_struct;
B2.m = 10;
B2.J = 0.31; % assumed (slender bar)
B2.r = [0.44212687; 0.60178858];
B2.p = 0;

B3 = Body_struct;
B3.m = 2;
B3.J = 0.0065; % assumed (slender bar)
B3.r = [0.55133463; 0.64329607];
B3.p = 0;

B4 = Body_struct;
B4.m = 2;
B4.J = 0.004; % assumed (slender bar)
B4.r = [0.61225784; 0.89109009];
B4.p = 0;

B5 = Body_struct;
B5.m = 8;
B5.J = 0.2; % assumed (slender bar)
B5.r = [1.53558308; 0.38018286];
B5.p = 0;

B6 = Body_struct;
B6.m = 10;
B6.J = 0.96; % assumed (m R^2)
B6.r = [0.34995625; 0.31];
B6.p = 0;

B7 = Body_struct;
B7.m = 10;
B7.J = 0.96; % assumed (m R^2)
B7.r = [1.8019937; 0.31];
B7.p = 0;

Bodies = [B1; B2; B3; B4; B5; B6; B7];
