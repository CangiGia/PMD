function ChainOpen(check, theta, theta_d, t)
% Recursive coordinate transformations
    include_global

if check <= 2
    % Update coordinates
    Coord_Float(1, theta(1:3));
    Coord_Rev(1, 2, theta(5));  % A1 -> A3
    Coord_Rev(5, 6, theta(4));  % D3 -> D2
    Coord_Rev(3, 4, theta(6));  % G1 -> G4
    Coord_Rev(9, 10, theta(7));  % I1 -> I5
    Coord_Rev(11, 12, theta(8));  % E2 -> E6
    Coord_Rev(13, 14, theta(9));  % K5 -> K7
    Update_P(7);
    Update_P(8);
    Update_P(15);
    Update_P(16);
    Update_P(17);
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
    Bmat(10:11, 1) = [1; 0];
    Bmat(10:11, 2) = [0; 1];
    Bmat(10:11, 3) = s_rot(Bodies(4).r - Bodies(1).r);
    Bmat(12, 3) = 1;
    Bmat(13:14, 1) = [1; 0];
    Bmat(13:14, 2) = [0; 1];
    Bmat(13:14, 3) = s_rot(Bodies(5).r - Bodies(1).r);
    Bmat(15, 3) = 1;
    Bmat(16:17, 1) = [1; 0];
    Bmat(16:17, 2) = [0; 1];
    Bmat(16:17, 3) = s_rot(Bodies(6).r - Bodies(1).r);
    Bmat(18, 3) = 1;
    Bmat(19:20, 1) = [1; 0];
    Bmat(19:20, 2) = [0; 1];
    Bmat(19:20, 3) = s_rot(Bodies(7).r - Bodies(1).r);
    Bmat(21, 3) = 1;
    Bmat(4:5, 5) = s_rot(Bodies(2).r - Points(2).rP);
    Bmat(6, 5) = 1;
    Bmat(7:8, 5) = s_rot(Bodies(3).r - Points(2).rP);
    Bmat(9, 5) = 1;
    Bmat(16:17, 5) = s_rot(Bodies(6).r - Points(2).rP);
    Bmat(18, 5) = 1;
    Bmat(4:5, 4) = s_rot(Bodies(2).r - Points(6).rP);
    Bmat(6, 4) = 1;
    Bmat(16:17, 4) = s_rot(Bodies(6).r - Points(6).rP);
    Bmat(18, 4) = 1;
    Bmat(10:11, 6) = s_rot(Bodies(4).r - Points(4).rP);
    Bmat(12, 6) = 1;
    Bmat(13:14, 7) = s_rot(Bodies(5).r - Points(10).rP);
    Bmat(15, 7) = 1;
    Bmat(19:20, 7) = s_rot(Bodies(7).r - Points(10).rP);
    Bmat(21, 7) = 1;
    Bmat(16:17, 8) = s_rot(Bodies(6).r - Points(12).rP);
    Bmat(18, 8) = 1;
    Bmat(19:20, 9) = s_rot(Bodies(7).r - Points(14).rP);
    Bmat(21, 9) = 1;

elseif check == 3
    % B_dot x theta_dot array
    Bd = zeros(nB3, 1);
    Bd(1:2) = s_rot(Bodies(1).r_d - Bodies(1).r_d)*theta_d(3);
    Bd(4:5) = s_rot(Bodies(2).r_d - Bodies(1).r_d)*theta_d(3) ...
        + s_rot(Bodies(2).r_d - Points(2).rP_d)*theta_d(5) ...
        + s_rot(Bodies(2).r_d - Points(6).rP_d)*theta_d(4);
    Bd(7:8) = s_rot(Bodies(3).r_d - Bodies(1).r_d)*theta_d(3) ...
        + s_rot(Bodies(3).r_d - Points(2).rP_d)*theta_d(5);
    Bd(10:11) = s_rot(Bodies(4).r_d - Bodies(1).r_d)*theta_d(3) ...
        + s_rot(Bodies(4).r_d - Points(4).rP_d)*theta_d(6);
    Bd(13:14) = s_rot(Bodies(5).r_d - Bodies(1).r_d)*theta_d(3) ...
        + s_rot(Bodies(5).r_d - Points(10).rP_d)*theta_d(7);
    Bd(16:17) = s_rot(Bodies(6).r_d - Bodies(1).r_d)*theta_d(3) ...
        + s_rot(Bodies(6).r_d - Points(2).rP_d)*theta_d(5) ...
        + s_rot(Bodies(6).r_d - Points(6).rP_d)*theta_d(4) ...
        + s_rot(Bodies(6).r_d - Points(12).rP_d)*theta_d(8);
    Bd(19:20) = s_rot(Bodies(7).r_d - Bodies(1).r_d)*theta_d(3) ...
        + s_rot(Bodies(7).r_d - Points(10).rP_d)*theta_d(7) ...
        + s_rot(Bodies(7).r_d - Points(14).rP_d)*theta_d(9);
end

end
