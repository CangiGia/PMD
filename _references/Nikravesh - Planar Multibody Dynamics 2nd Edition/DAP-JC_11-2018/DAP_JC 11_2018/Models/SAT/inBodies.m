function inBodies
    include_global

% Satellite solar panels, phase 1 (Nikravesh, Section 15.17): joint-coordinate model

B1 = Body_struct;
B1.m = 1;
B1.J = 0.03; % assumed (thin plate)
B1.theta = 2.0943951;

B2 = Body_struct;
B2.m = 0.035;
B2.J = 0.00035729; % assumed (slender bar)
B2.theta = -2.11955084;

B3 = Body_struct;
B3.m = 0.035;
B3.J = 0.00035729; % assumed (slender bar)
B3.theta = 2.14470658;

B4 = Body_struct;
B4.m = 0.03;
B4.J = 0.000225; % assumed (slender bar)
B4.theta = -0.04025853;

B5 = Body_struct;
B5.m = 0.03;
B5.J = 0.000225; % assumed (slender bar)
B5.theta = -2.57779055;

Bodies = [B1; B2; B3; B4; B5];
