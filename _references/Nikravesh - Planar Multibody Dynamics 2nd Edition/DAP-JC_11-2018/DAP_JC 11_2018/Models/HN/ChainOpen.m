function ChainOpen(check, theta, theta_d, t)
% Recursive coordinate transformations
    include_global

if check <= 2
    % Update coordinates
    Coord_Tran(1, 1, 2, 2, theta(1));
    Coord_Rev(3, 4, theta(2));  % R0a -> R0b
    Coord_Rev(5, 6, theta(3));  % R1a -> R1b
    Coord_Rev(7, 8, theta(4));  % R2a -> R2b
    Coord_Rev(9, 10, theta(5));  % R3a -> R3b
    Coord_Rev(11, 12, theta(6));  % R4a -> R4b
    Coord_Rev(13, 14, theta(7));  % R5a -> R5b
    Coord_Rev(15, 16, theta(8));  % R6a -> R6b
    Coord_Rev(17, 18, theta(9));  % R7a -> R7b
    Coord_Float(10, theta(10:12));
    Update_P(19);
    Update_P(20);
    Update_P(21);
    Update_P(22);
    Update_P(23);
end

if check == 2
    % B matrix
    Bmat = zeros(nB3, nJC);
    Bmat(1:2, 1) = Uvectors(2).u;
    Bmat(4:5, 1) = Uvectors(2).u;
    Bmat(7:8, 1) = Uvectors(2).u;
    Bmat(10:11, 1) = Uvectors(2).u;
    Bmat(13:14, 1) = Uvectors(2).u;
    Bmat(16:17, 1) = Uvectors(2).u;
    Bmat(19:20, 1) = Uvectors(2).u;
    Bmat(22:23, 1) = Uvectors(2).u;
    Bmat(25:26, 1) = Uvectors(2).u;
    Bmat(4:5, 2) = s_rot(Bodies(2).r - Points(4).rP);
    Bmat(6, 2) = 1;
    Bmat(7:8, 2) = s_rot(Bodies(3).r - Points(4).rP);
    Bmat(9, 2) = 1;
    Bmat(10:11, 2) = s_rot(Bodies(4).r - Points(4).rP);
    Bmat(12, 2) = 1;
    Bmat(13:14, 2) = s_rot(Bodies(5).r - Points(4).rP);
    Bmat(15, 2) = 1;
    Bmat(16:17, 2) = s_rot(Bodies(6).r - Points(4).rP);
    Bmat(18, 2) = 1;
    Bmat(19:20, 2) = s_rot(Bodies(7).r - Points(4).rP);
    Bmat(21, 2) = 1;
    Bmat(22:23, 2) = s_rot(Bodies(8).r - Points(4).rP);
    Bmat(24, 2) = 1;
    Bmat(25:26, 2) = s_rot(Bodies(9).r - Points(4).rP);
    Bmat(27, 2) = 1;
    Bmat(7:8, 3) = s_rot(Bodies(3).r - Points(6).rP);
    Bmat(9, 3) = 1;
    Bmat(10:11, 3) = s_rot(Bodies(4).r - Points(6).rP);
    Bmat(12, 3) = 1;
    Bmat(13:14, 3) = s_rot(Bodies(5).r - Points(6).rP);
    Bmat(15, 3) = 1;
    Bmat(16:17, 3) = s_rot(Bodies(6).r - Points(6).rP);
    Bmat(18, 3) = 1;
    Bmat(19:20, 3) = s_rot(Bodies(7).r - Points(6).rP);
    Bmat(21, 3) = 1;
    Bmat(22:23, 3) = s_rot(Bodies(8).r - Points(6).rP);
    Bmat(24, 3) = 1;
    Bmat(25:26, 3) = s_rot(Bodies(9).r - Points(6).rP);
    Bmat(27, 3) = 1;
    Bmat(10:11, 4) = s_rot(Bodies(4).r - Points(8).rP);
    Bmat(12, 4) = 1;
    Bmat(13:14, 4) = s_rot(Bodies(5).r - Points(8).rP);
    Bmat(15, 4) = 1;
    Bmat(16:17, 4) = s_rot(Bodies(6).r - Points(8).rP);
    Bmat(18, 4) = 1;
    Bmat(19:20, 4) = s_rot(Bodies(7).r - Points(8).rP);
    Bmat(21, 4) = 1;
    Bmat(22:23, 4) = s_rot(Bodies(8).r - Points(8).rP);
    Bmat(24, 4) = 1;
    Bmat(25:26, 4) = s_rot(Bodies(9).r - Points(8).rP);
    Bmat(27, 4) = 1;
    Bmat(13:14, 5) = s_rot(Bodies(5).r - Points(10).rP);
    Bmat(15, 5) = 1;
    Bmat(16:17, 5) = s_rot(Bodies(6).r - Points(10).rP);
    Bmat(18, 5) = 1;
    Bmat(19:20, 5) = s_rot(Bodies(7).r - Points(10).rP);
    Bmat(21, 5) = 1;
    Bmat(22:23, 5) = s_rot(Bodies(8).r - Points(10).rP);
    Bmat(24, 5) = 1;
    Bmat(25:26, 5) = s_rot(Bodies(9).r - Points(10).rP);
    Bmat(27, 5) = 1;
    Bmat(16:17, 6) = s_rot(Bodies(6).r - Points(12).rP);
    Bmat(18, 6) = 1;
    Bmat(19:20, 6) = s_rot(Bodies(7).r - Points(12).rP);
    Bmat(21, 6) = 1;
    Bmat(22:23, 6) = s_rot(Bodies(8).r - Points(12).rP);
    Bmat(24, 6) = 1;
    Bmat(25:26, 6) = s_rot(Bodies(9).r - Points(12).rP);
    Bmat(27, 6) = 1;
    Bmat(19:20, 7) = s_rot(Bodies(7).r - Points(14).rP);
    Bmat(21, 7) = 1;
    Bmat(22:23, 7) = s_rot(Bodies(8).r - Points(14).rP);
    Bmat(24, 7) = 1;
    Bmat(25:26, 7) = s_rot(Bodies(9).r - Points(14).rP);
    Bmat(27, 7) = 1;
    Bmat(22:23, 8) = s_rot(Bodies(8).r - Points(16).rP);
    Bmat(24, 8) = 1;
    Bmat(25:26, 8) = s_rot(Bodies(9).r - Points(16).rP);
    Bmat(27, 8) = 1;
    Bmat(25:26, 9) = s_rot(Bodies(9).r - Points(18).rP);
    Bmat(27, 9) = 1;
    Bmat(28:29, 10) = [1; 0];
    Bmat(28:29, 11) = [0; 1];
    Bmat(28:29, 12) = s_rot(Bodies(10).r - Bodies(10).r);
    Bmat(30, 12) = 1;

elseif check == 3
    % B_dot x theta_dot array
    Bd = zeros(nB3, 1);
    Bd(4:5) = s_rot(Bodies(2).r_d - Points(4).rP_d)*theta_d(2);
    Bd(7:8) = s_rot(Bodies(3).r_d - Points(4).rP_d)*theta_d(2) ...
        + s_rot(Bodies(3).r_d - Points(6).rP_d)*theta_d(3);
    Bd(10:11) = s_rot(Bodies(4).r_d - Points(4).rP_d)*theta_d(2) ...
        + s_rot(Bodies(4).r_d - Points(6).rP_d)*theta_d(3) ...
        + s_rot(Bodies(4).r_d - Points(8).rP_d)*theta_d(4);
    Bd(13:14) = s_rot(Bodies(5).r_d - Points(4).rP_d)*theta_d(2) ...
        + s_rot(Bodies(5).r_d - Points(6).rP_d)*theta_d(3) ...
        + s_rot(Bodies(5).r_d - Points(8).rP_d)*theta_d(4) ...
        + s_rot(Bodies(5).r_d - Points(10).rP_d)*theta_d(5);
    Bd(16:17) = s_rot(Bodies(6).r_d - Points(4).rP_d)*theta_d(2) ...
        + s_rot(Bodies(6).r_d - Points(6).rP_d)*theta_d(3) ...
        + s_rot(Bodies(6).r_d - Points(8).rP_d)*theta_d(4) ...
        + s_rot(Bodies(6).r_d - Points(10).rP_d)*theta_d(5) ...
        + s_rot(Bodies(6).r_d - Points(12).rP_d)*theta_d(6);
    Bd(19:20) = s_rot(Bodies(7).r_d - Points(4).rP_d)*theta_d(2) ...
        + s_rot(Bodies(7).r_d - Points(6).rP_d)*theta_d(3) ...
        + s_rot(Bodies(7).r_d - Points(8).rP_d)*theta_d(4) ...
        + s_rot(Bodies(7).r_d - Points(10).rP_d)*theta_d(5) ...
        + s_rot(Bodies(7).r_d - Points(12).rP_d)*theta_d(6) ...
        + s_rot(Bodies(7).r_d - Points(14).rP_d)*theta_d(7);
    Bd(22:23) = s_rot(Bodies(8).r_d - Points(4).rP_d)*theta_d(2) ...
        + s_rot(Bodies(8).r_d - Points(6).rP_d)*theta_d(3) ...
        + s_rot(Bodies(8).r_d - Points(8).rP_d)*theta_d(4) ...
        + s_rot(Bodies(8).r_d - Points(10).rP_d)*theta_d(5) ...
        + s_rot(Bodies(8).r_d - Points(12).rP_d)*theta_d(6) ...
        + s_rot(Bodies(8).r_d - Points(14).rP_d)*theta_d(7) ...
        + s_rot(Bodies(8).r_d - Points(16).rP_d)*theta_d(8);
    Bd(25:26) = s_rot(Bodies(9).r_d - Points(4).rP_d)*theta_d(2) ...
        + s_rot(Bodies(9).r_d - Points(6).rP_d)*theta_d(3) ...
        + s_rot(Bodies(9).r_d - Points(8).rP_d)*theta_d(4) ...
        + s_rot(Bodies(9).r_d - Points(10).rP_d)*theta_d(5) ...
        + s_rot(Bodies(9).r_d - Points(12).rP_d)*theta_d(6) ...
        + s_rot(Bodies(9).r_d - Points(14).rP_d)*theta_d(7) ...
        + s_rot(Bodies(9).r_d - Points(16).rP_d)*theta_d(8) ...
        + s_rot(Bodies(9).r_d - Points(18).rP_d)*theta_d(9);
    Bd(28:29) = s_rot(Bodies(10).r_d - Bodies(10).r_d)*theta_d(12);
end

end
