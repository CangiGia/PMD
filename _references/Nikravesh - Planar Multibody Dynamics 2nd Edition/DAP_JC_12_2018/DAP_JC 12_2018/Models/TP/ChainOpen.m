function ChainOpen(check, theta, theta_d, t)
% Recursive coordinate transformations
    include_global
    
if check <= 2
    % Update coordinates
    Coord_Rev(1, 2, theta(1));
    Coord_Rev(3, 5, theta(2));
    Coord_Rev(4, 6, theta(3));
end
    
if check == 2
    % B matrix
    d11 =  Bodies(1).r;
    d21 =  Bodies(2).r;
    d31 =  Bodies(3).r;
    d22 = -Points(5).sP;
    d33 = -Points(6).sP;
    z2 = [0; 0];
    Bmat = [s_rot(d11)  z2    z2
            1    0    0
            s_rot(d21)  s_rot(d22)    z2
            1    1    0
            s_rot(d31)  z2   s_rot(d33)
            1    0    1];
            
elseif check == 3
    % B_dot x theta_dot array
    
    d11d =  Bodies(1).r_d;
    d21d =  Bodies(2).r_d;
    d31d =  Bodies(3).r_d;
    d22d = -Points(5).sP_d;
    d33d = -Points(6).sP_d;
    
    Bd = [s_rot(d11d)*theta_d(1)
            0
            s_rot(d21d)*theta_d(1) + s_rot(d22d)*theta_d(2)
            0
            s_rot(d31d)*theta_d(1) + s_rot(d33d)*theta_d(3)
            0];
 end

end
    
     