function ChainOpen(check, theta, theta_d, t)
% Recursive coordinate transformations
    include_global
    
if check <= 2
    % Update coordinates
    Coord_Float(1, theta(1:3));
    Update_P(1);
    Update_P(2);
end
    
if check == 2
    % B matrix
    
    Bmat = [eye(3)];
            
elseif check == 3
    % B_dot x theta_dot array
    
    Bd = [0
          0
          0];
end

end
    
     