function ChainOpen(check, theta, theta_d, t)
% Recursive coordinate transformations
    include_global

if check <= 2
    % Update coordinates
    Coord_Rev(1, 3, theta(1));  % A0 -> A1
    Coord_Rev(4, 5, theta(2));  % C1 -> C2
    Coord_Rev(2, 10, theta(3));  % B0 -> B3
    Coord_Rev(9, 11, theta(4));  % E3 -> E4
    Coord_Rev(12, 13, theta(5));  % F4 -> F5
    Update_P(6);
    Update_P(7);
    Update_P(8);
    Update_P(14);
end

if check == 2
    % B matrix
    Bmat = zeros(nB3, nJC);
    Bmat(1:2, 1) = s_rot(Bodies(1).r - Points(3).rP);
    Bmat(3, 1) = 1;
    Bmat(4:5, 1) = s_rot(Bodies(2).r - Points(3).rP);
    Bmat(6, 1) = 1;
    Bmat(4:5, 2) = s_rot(Bodies(2).r - Points(5).rP);
    Bmat(6, 2) = 1;
    Bmat(7:8, 3) = s_rot(Bodies(3).r - Points(10).rP);
    Bmat(9, 3) = 1;
    Bmat(10:11, 3) = s_rot(Bodies(4).r - Points(10).rP);
    Bmat(12, 3) = 1;
    Bmat(13:14, 3) = s_rot(Bodies(5).r - Points(10).rP);
    Bmat(15, 3) = 1;
    Bmat(10:11, 4) = s_rot(Bodies(4).r - Points(11).rP);
    Bmat(12, 4) = 1;
    Bmat(13:14, 4) = s_rot(Bodies(5).r - Points(11).rP);
    Bmat(15, 4) = 1;
    Bmat(13:14, 5) = s_rot(Bodies(5).r - Points(13).rP);
    Bmat(15, 5) = 1;

elseif check == 3
    % B_dot x theta_dot array
    Bd = zeros(nB3, 1);
    Bd(1:2) = s_rot(Bodies(1).r_d - Points(3).rP_d)*theta_d(1);
    Bd(4:5) = s_rot(Bodies(2).r_d - Points(3).rP_d)*theta_d(1) ...
        + s_rot(Bodies(2).r_d - Points(5).rP_d)*theta_d(2);
    Bd(7:8) = s_rot(Bodies(3).r_d - Points(10).rP_d)*theta_d(3);
    Bd(10:11) = s_rot(Bodies(4).r_d - Points(10).rP_d)*theta_d(3) ...
        + s_rot(Bodies(4).r_d - Points(11).rP_d)*theta_d(4);
    Bd(13:14) = s_rot(Bodies(5).r_d - Points(10).rP_d)*theta_d(3) ...
        + s_rot(Bodies(5).r_d - Points(11).rP_d)*theta_d(4) ...
        + s_rot(Bodies(5).r_d - Points(13).rP_d)*theta_d(5);
end

end
