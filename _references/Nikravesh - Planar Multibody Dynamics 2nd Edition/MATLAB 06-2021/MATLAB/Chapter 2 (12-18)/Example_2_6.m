% Example 2.6
% Define square matrix D
    D = [-1  -1.64   2.14
          0   2.51  -0.51
          1   0      0   ]
% Define right-hand-side array
    rhs = [0  0  2*pi]'
% Solve for the unknowns using inverse of D
    ans_1 = inv(D)*rhs
% Solve for the unknowns using backslash
    ans_2 = D\rhs
