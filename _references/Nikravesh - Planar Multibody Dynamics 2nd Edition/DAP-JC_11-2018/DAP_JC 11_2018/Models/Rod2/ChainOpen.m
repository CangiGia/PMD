function ChainOpen(check, theta, theta_d, t)
% Recursive coordinate transformations
    include_global
    
if check <= 2
    % Update coordinates
    Coord_Float(1, theta(1:3));
    Coord_Rev(2, 3, theta(4));
    Update_P(1);
    Update_P(4);
end
    
if check == 2
    % B matrix
    d21 =  Bodies(2).r - Bodies(1).r;
    d22 = -Points(3).sP;

    Bmat = [eye(3) [0; 0; 0]
            eye(2) s_rot(d21)  s_rot(d22)
            0  0   1           1];
            
elseif check == 3
    % B_dot x theta_dot array
    d21d =  Bodies(2).r_d - Bodies(1).r_d;
    d22d = -Points(3).sP_d;

    Bd = [0
          0
          0
          s_rot(d21d*theta_d(3) + d22d*theta_d(4))
          0];
end

end
    
     