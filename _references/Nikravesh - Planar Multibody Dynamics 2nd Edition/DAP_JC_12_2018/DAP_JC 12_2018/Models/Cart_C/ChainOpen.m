function ChainOpen(check, theta, theta_d, t)
% Recursive coordinate transformations
    include_global
    
if check <= 2
    % Update coordinates
    R = 0.1; a = 0.6;
    Bodies(1).r = [-R*theta(1); R];
    Bodies(1).p = theta(1);
    Bodies(1).A = Matrix_A(theta(1));
    Update_P(1);
    Bodies(2).p = theta(2);
    Bodies(2).A = Matrix_A(theta(2));
    Points(2).sP = Bodies(2).A*Points(2).sPlocal;
    Points(2).sP_r = s_rot(Points(2).sP);
    Points(2).rP = Points(1).rP;
    Bodies(2).r = Points(2).rP - Points(2).sP;
    Update_P(3);
    Bodies(3).r = [a-R*theta(3); R];
    Bodies(3).p = theta(3);
    Bodies(3).A = Matrix_A(theta(3));
    Update_P(4);
end
    
if check == 2
    % B matrix
    Rux = [0.1; 0];
    z2 = [0; 0];
    Bmat = [-Rux   z2    z2
             1     0     0
            -Rux  -Points(2).sP_r    z2
             0     1     0
             z2    z2   -Rux
             0     0     1];
           
elseif check == 3
    % B_dot x theta_dot array
    z2 = [0; 0];
    Bd = [z2
          0
          Points(2).sP_d*theta_d(2)
          0
          z2
          0];
 end

end
    
     