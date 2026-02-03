% M-file Gees

%  Unit vectors for mass centers 
    [u, ur]   = uVectors(theta);
    [uG, urG] = uVectors(theta + gamma);
% Coordinates of G1, G2, G3
    RG1(i,:) = LG(1)*uG{1}'; 
    RG2(i,:) = L(1)*u{1}' + LG(2)*uG{2}'; 
    RG3(i,:) = RQ' + LG(3)*uG{3}';
% Velocity of G1, G2, G3
    VG1(i,:) = LG(1)*urG{1}'*omegas(i,1); 
    VG2(i,:) = L(1)*ur{1}'*omegas(i,1) + LG(2)*urG{2}'*omegas(i,2); 
    VG3(i,:) = LG(3)*urG{3}'*omegas(i,3);
% Acceleration of G1, G2, G3
    AG1(i,:) = - LG(1)*uG{1}'*omegas(i,1)^2; 
    AG2(i,:) = AA(i,:) + LG(2)*(urG{2}'*alphas(i,2) - uG{2}'*omegas(i,2)^2);
    AG3(i,:) = LG(3)*(urG{3}'*alphas(i,3) - uG{3}'*omegas(i,3)^2);
