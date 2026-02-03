% M-file Pees

% Unit vectors
    [uP, urP] = uVectors(theta + beta);
% Coordinates of A, B, C and P
    RA(i,:)  = L(1)*u{1}'; 
    RB(i,:)  = RA(i,:) + L(2)*u{2}';
    RP1(i,:) = LP(1)*uP{1}';
    RP2(i,:) = RA(i,:) + LP(2)*uP{2}';
    RP3(i,:) = RQ' + LP(3)*uP{3}';
% Velocity of A, B, C and P
    VA(i,:)  = L(1)*ur{1}'*omega(1); 
    VB(i,:)  = L(3)*ur{3}'*omega(3); 
    VP1(i,:) = LP(1)*urP{1}'*omega(1);
    VP2(i,:) = VA(i,:) + LP(2)*urP{2}'*omega(2);
    VP3(i,:) = LP(3)*urP{3}'*omega(3);
% Acceleration of A, B, C and P
    AA(i,:)  = -RA(i,:)'*omega(1)^2; 
    AB(i,:)  = L(3)*(ur{3}'*alpha(3) - u{3}'*omega(3)^2); 
    AP1(i,:) = -LP(1)*uP{1}'*omega(1)^2;
    AP2(i,:) = AA(i,:) + LP(2)*(urP{2}'*alpha(2) - uP{2}'*omega(2)^2);
    AP3(i,:) = LP(3)*(urP{3}'*alpha(3) - uP{3}'*omega(3)^2);
