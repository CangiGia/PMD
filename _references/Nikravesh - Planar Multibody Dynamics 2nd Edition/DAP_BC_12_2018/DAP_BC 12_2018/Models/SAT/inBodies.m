function inBodies
    include_global

% Satellite solar panels, phase 1 (Nikravesh, Section 15.17): body-coordinate model

B1 = Body_struct;
B1.m = 1;
B1.J = 0.03; % assumed (thin plate)
B1.r = [-0.48464102; -0.31980762];
B1.p = 2.0943951;

B2 = Body_struct;
B2.m = 0.035;
B2.J = 0.00035729; % assumed (slender bar)
B2.r = [-0.29237588; -0.31420941];
B2.p = -0.02515574;

B3 = Body_struct;
B3.m = 0.035;
B3.J = 0.00035729; % assumed (slender bar)
B3.r = [-0.20871562; -0.1693056];
B3.p = 2.11955084;

B4 = Body_struct;
B4.m = 0.03;
B4.J = 0.000225; % assumed (slender bar)
B4.r = [-0.19046086; -0.18758954];
B4.p = 2.07929231;

B5 = Body_struct;
B5.m = 0.03;
B5.J = 0.000225; % assumed (slender bar)
B5.r = [-0.13174523; -0.12828394];
B5.p = -0.49849823;

Bodies = [B1; B2; B3; B4; B5];
