function ChainOpen(check, theta, theta_d, t)
% Recursive coordinate transformations
    include_global
    
if check <= 2
    % Update coordinates
    Coord_Rev(7, 1, theta(1));
    Coord_Rev(2, 3, theta(2));
    Update_P(5);
    Update_U(1);
end
    
if check == 2
    % B matrix
    d11 = -Points(1).sP;
    d21 =  Bodies(2).r - Points(1).rP;
    d22 = -Points(3).sP;
    z2 = [0; 0];
    Bmat = [s_rot(d11)  z2   
            1  0
            s_rot(d21)  s_rot(d22) 
            1           1];
            
elseif check == 3
    % B_dot x theta_dot array
    d11d =  Bodies(1).r_d;
    d21d =  Bodies(2).r_d;
    d22d = -Points(3).sP_d;

    Bd = [s_rot(d11d)*theta_d(1)
          0
          s_rot(d21d*theta_d(1) + d22d*theta_d(2))
          0];
end

end
    
     