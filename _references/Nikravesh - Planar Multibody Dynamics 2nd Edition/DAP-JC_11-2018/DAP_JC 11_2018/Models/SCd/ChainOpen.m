function ChainOpen(check, theta, theta_d, t)
% Recursive coordinate transformations
    include_global
    
if check <= 2
    % Update coordinates
    Coord_Rev(1, 2, theta(1));
    Coord_Rev(3, 4, theta(2));
    Coord_Tran(1, 1, 6, 2, theta(3));
    Update_P(5);
end
    
if check == 2
    % B matrix
    d11 = -Points(2).sP;
    d21 =  Bodies(2).r;
    d22 = -Points(4).sP;
    u3 = Uvectors(2).u;
    z2 = [0; 0];
    Bmat = [s_rot(d11)  z2          z2
            1           0           0
            s_rot(d21)  s_rot(d22)  z2
            1           1           0
            z2          z2          u3
            0           0           0];
            
elseif check == 3
    % B_dot x theta_dot array
    d11d = -Points(2).sP_d;
    d21d =  Bodies(2).r_d;
    d22d = -Points(4).sP_d;

    Bd = [s_rot(d11d)*theta_d(1)
          0
          s_rot(d21d*theta_d(1) + d22d*theta_d(2))
          0
          0
          0
          0];
end

end
    
     