% M-file moms

% Moment arms w.r.t the mass centers 
    sO1 =          - RG1(i,:); sO1r = s_rot(sO1);
    sA1 = RA(i,:)  - RG1(i,:); sA1r = s_rot(sA1);
    sA2 = RA(i,:)  - RG2(i,:); sA2r = s_rot(sA2);
    sB2 = RB(i,:)  - RG2(i,:); sB2r = s_rot(sB2);
    sB3 = RB(i,:)  - RG3(i,:); sB3r = s_rot(sB3);
    sQ3 = RQ'      - RG3(i,:); sQ3r = s_rot(sQ3);
    sP1 = RP1(i,:) - RG1(i,:); sP1r = s_rot(sP1);
    sP2 = RP2(i,:) - RG2(i,:); sP2r = s_rot(sP2);
    sP3 = RP3(i,:) - RG3(i,:); sP3r = s_rot(sP3);

% Jacobian
    I2 = eye(2); Z23 = zeros(2,3);
    O1 = [ I2  sO1r]; A1 = [-I2 -sA1r];
    A2 = [ I2  sA2r]; B2 = [-I2 -sB2r];
    B3 = [ I2  sB3r]; Q3 = [ I2  sQ3r];
    T1 = [0 0 1]; Z13 = zeros(1,3);
    D  = [O1  Z23 Z23
          A1  A2  Z23
          Z23 B2  B3
          Z23 Z23 Q3
          T1  Z13 Z13];
