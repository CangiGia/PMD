function ChainOpen(check, theta, theta_d, t)
% Recursive coordinate transformations
    include_global
    
% -----------------------------------------    
if check <= 2
    % Update coordinates
    Coord_Tran(1, 1, 2, 2, theta(1));
    Coord_Rev(2, 3, theta(2));
end
    
% -----------------------------------------    
if check == 2
    % B matrix
    u1 = Uvectors(2).u;
    d22 = -Points(3).sP;
    z2 = [0; 0];
    Bmat = [u1  z2
            0   0
            u1  s_rot(d22)
            0   1];
            
% -----------------------------------------    
elseif check == 3
    % B_dot x theta_dot array
    
    d22d = -Points(3).sP_d;

    Bd = [0
          0
          0
          s_rot(d22d)*theta_d(2)
          0];
end

end
    
     