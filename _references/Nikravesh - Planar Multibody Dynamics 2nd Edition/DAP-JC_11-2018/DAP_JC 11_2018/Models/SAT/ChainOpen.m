function ChainOpen(check, theta, theta_d, t)
% Recursive coordinate transformations
    include_global

if check <= 2
    % Update coordinates
    Coord_Rev(1, 4, theta(1));  % B0 -> B1
    Coord_Rev(5, 6, theta(2));  % A1 -> A2
    Coord_Rev(7, 8, theta(3));  % S12 -> S13
    Coord_Rev(9, 11, theta(4));  % S13b -> S14
    Coord_Rev(12, 13, theta(5));  % L14 -> L15
    Update_P(10);
    Update_P(14);
end

if check == 2
    % B matrix
    Bmat = zeros(nB3, nJC);
    Bmat(1:2, 1) = s_rot(Bodies(1).r - Points(4).rP);
    Bmat(3, 1) = 1;
    Bmat(4:5, 1) = s_rot(Bodies(2).r - Points(4).rP);
    Bmat(6, 1) = 1;
    Bmat(7:8, 1) = s_rot(Bodies(3).r - Points(4).rP);
    Bmat(9, 1) = 1;
    Bmat(10:11, 1) = s_rot(Bodies(4).r - Points(4).rP);
    Bmat(12, 1) = 1;
    Bmat(13:14, 1) = s_rot(Bodies(5).r - Points(4).rP);
    Bmat(15, 1) = 1;
    Bmat(4:5, 2) = s_rot(Bodies(2).r - Points(6).rP);
    Bmat(6, 2) = 1;
    Bmat(7:8, 2) = s_rot(Bodies(3).r - Points(6).rP);
    Bmat(9, 2) = 1;
    Bmat(10:11, 2) = s_rot(Bodies(4).r - Points(6).rP);
    Bmat(12, 2) = 1;
    Bmat(13:14, 2) = s_rot(Bodies(5).r - Points(6).rP);
    Bmat(15, 2) = 1;
    Bmat(7:8, 3) = s_rot(Bodies(3).r - Points(8).rP);
    Bmat(9, 3) = 1;
    Bmat(10:11, 3) = s_rot(Bodies(4).r - Points(8).rP);
    Bmat(12, 3) = 1;
    Bmat(13:14, 3) = s_rot(Bodies(5).r - Points(8).rP);
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
    Bd(1:2) = s_rot(Bodies(1).r_d - Points(4).rP_d)*theta_d(1);
    Bd(4:5) = s_rot(Bodies(2).r_d - Points(4).rP_d)*theta_d(1) ...
        + s_rot(Bodies(2).r_d - Points(6).rP_d)*theta_d(2);
    Bd(7:8) = s_rot(Bodies(3).r_d - Points(4).rP_d)*theta_d(1) ...
        + s_rot(Bodies(3).r_d - Points(6).rP_d)*theta_d(2) ...
        + s_rot(Bodies(3).r_d - Points(8).rP_d)*theta_d(3);
    Bd(10:11) = s_rot(Bodies(4).r_d - Points(4).rP_d)*theta_d(1) ...
        + s_rot(Bodies(4).r_d - Points(6).rP_d)*theta_d(2) ...
        + s_rot(Bodies(4).r_d - Points(8).rP_d)*theta_d(3) ...
        + s_rot(Bodies(4).r_d - Points(11).rP_d)*theta_d(4);
    Bd(13:14) = s_rot(Bodies(5).r_d - Points(4).rP_d)*theta_d(1) ...
        + s_rot(Bodies(5).r_d - Points(6).rP_d)*theta_d(2) ...
        + s_rot(Bodies(5).r_d - Points(8).rP_d)*theta_d(3) ...
        + s_rot(Bodies(5).r_d - Points(11).rP_d)*theta_d(4) ...
        + s_rot(Bodies(5).r_d - Points(13).rP_d)*theta_d(5);
end

end
