function ChainOpen(check, theta, theta_d, t)
% Recursive coordinate transformations
    include_global

if check <= 2
    % Update coordinates
    Coord_Float(1, theta(1:3));
    Coord_Rev(1, 2, theta(4));  % O1 -> O2
    Coord_Rev(3, 4, theta(5));  % A2 -> A3
    Coord_Float(4, theta(6:8));
    Update_P(5);
    Update_P(6);
    Update_P(7);
    Update_P(8);
    Update_U(1);
end

if check == 2
    % B matrix
    Bmat = zeros(nB3, nJC);
    Bmat(1:2, 1) = [1; 0];
    Bmat(1:2, 2) = [0; 1];
    Bmat(1:2, 3) = s_rot(Bodies(1).r - Bodies(1).r);
    Bmat(3, 3) = 1;
    Bmat(4:5, 1) = [1; 0];
    Bmat(4:5, 2) = [0; 1];
    Bmat(4:5, 3) = s_rot(Bodies(2).r - Bodies(1).r);
    Bmat(6, 3) = 1;
    Bmat(7:8, 1) = [1; 0];
    Bmat(7:8, 2) = [0; 1];
    Bmat(7:8, 3) = s_rot(Bodies(3).r - Bodies(1).r);
    Bmat(9, 3) = 1;
    Bmat(4:5, 4) = s_rot(Bodies(2).r - Points(2).rP);
    Bmat(6, 4) = 1;
    Bmat(7:8, 4) = s_rot(Bodies(3).r - Points(2).rP);
    Bmat(9, 4) = 1;
    Bmat(7:8, 5) = s_rot(Bodies(3).r - Points(4).rP);
    Bmat(9, 5) = 1;
    Bmat(10:11, 6) = [1; 0];
    Bmat(10:11, 7) = [0; 1];
    Bmat(10:11, 8) = s_rot(Bodies(4).r - Bodies(4).r);
    Bmat(12, 8) = 1;

elseif check == 3
    % B_dot x theta_dot array
    Bd = zeros(nB3, 1);
    Bd(1:2) = s_rot(Bodies(1).r_d - Bodies(1).r_d)*theta_d(3);
    Bd(4:5) = s_rot(Bodies(2).r_d - Bodies(1).r_d)*theta_d(3) ...
        + s_rot(Bodies(2).r_d - Points(2).rP_d)*theta_d(4);
    Bd(7:8) = s_rot(Bodies(3).r_d - Bodies(1).r_d)*theta_d(3) ...
        + s_rot(Bodies(3).r_d - Points(2).rP_d)*theta_d(4) ...
        + s_rot(Bodies(3).r_d - Points(4).rP_d)*theta_d(5);
    Bd(10:11) = s_rot(Bodies(4).r_d - Bodies(4).r_d)*theta_d(8);
end

end
