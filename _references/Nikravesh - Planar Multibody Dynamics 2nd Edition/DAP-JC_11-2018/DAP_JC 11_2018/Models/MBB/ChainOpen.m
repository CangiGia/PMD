function ChainOpen(check, theta, theta_d, t)
% Recursive coordinate transformations
    include_global

if check <= 2
    % Update coordinates
    Coord_Float(1, theta(1:3));
    Coord_Rev(1, 2, theta(4));  % O1 -> O2
    Coord_Rev(3, 4, theta(5));  % A2 -> A3
    Coord_Rev(7, 8, theta(6));  % B3 -> B4
    Coord_Rev(5, 6, theta(7));  % A2w -> A5
    Coord_Float(6, theta(8:10));
    Update_P(9);
    Update_P(10);
    Update_P(11);
    Update_P(12);
    Update_P(13);
    Update_P(14);
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
    Bmat(10:11, 1) = [1; 0];
    Bmat(10:11, 2) = [0; 1];
    Bmat(10:11, 3) = s_rot(Bodies(4).r - Bodies(1).r);
    Bmat(12, 3) = 1;
    Bmat(13:14, 1) = [1; 0];
    Bmat(13:14, 2) = [0; 1];
    Bmat(13:14, 3) = s_rot(Bodies(5).r - Bodies(1).r);
    Bmat(15, 3) = 1;
    Bmat(4:5, 4) = s_rot(Bodies(2).r - Points(2).rP);
    Bmat(6, 4) = 1;
    Bmat(7:8, 4) = s_rot(Bodies(3).r - Points(2).rP);
    Bmat(9, 4) = 1;
    Bmat(10:11, 4) = s_rot(Bodies(4).r - Points(2).rP);
    Bmat(12, 4) = 1;
    Bmat(13:14, 4) = s_rot(Bodies(5).r - Points(2).rP);
    Bmat(15, 4) = 1;
    Bmat(7:8, 5) = s_rot(Bodies(3).r - Points(4).rP);
    Bmat(9, 5) = 1;
    Bmat(10:11, 5) = s_rot(Bodies(4).r - Points(4).rP);
    Bmat(12, 5) = 1;
    Bmat(10:11, 6) = s_rot(Bodies(4).r - Points(8).rP);
    Bmat(12, 6) = 1;
    Bmat(13:14, 7) = s_rot(Bodies(5).r - Points(6).rP);
    Bmat(15, 7) = 1;
    Bmat(16:17, 8) = [1; 0];
    Bmat(16:17, 9) = [0; 1];
    Bmat(16:17, 10) = s_rot(Bodies(6).r - Bodies(6).r);
    Bmat(18, 10) = 1;

elseif check == 3
    % B_dot x theta_dot array
    Bd = zeros(nB3, 1);
    Bd(1:2) = s_rot(Bodies(1).r_d - Bodies(1).r_d)*theta_d(3);
    Bd(4:5) = s_rot(Bodies(2).r_d - Bodies(1).r_d)*theta_d(3) ...
        + s_rot(Bodies(2).r_d - Points(2).rP_d)*theta_d(4);
    Bd(7:8) = s_rot(Bodies(3).r_d - Bodies(1).r_d)*theta_d(3) ...
        + s_rot(Bodies(3).r_d - Points(2).rP_d)*theta_d(4) ...
        + s_rot(Bodies(3).r_d - Points(4).rP_d)*theta_d(5);
    Bd(10:11) = s_rot(Bodies(4).r_d - Bodies(1).r_d)*theta_d(3) ...
        + s_rot(Bodies(4).r_d - Points(2).rP_d)*theta_d(4) ...
        + s_rot(Bodies(4).r_d - Points(4).rP_d)*theta_d(5) ...
        + s_rot(Bodies(4).r_d - Points(8).rP_d)*theta_d(6);
    Bd(13:14) = s_rot(Bodies(5).r_d - Bodies(1).r_d)*theta_d(3) ...
        + s_rot(Bodies(5).r_d - Points(2).rP_d)*theta_d(4) ...
        + s_rot(Bodies(5).r_d - Points(6).rP_d)*theta_d(7);
    Bd(16:17) = s_rot(Bodies(6).r_d - Bodies(6).r_d)*theta_d(10);
end

end
