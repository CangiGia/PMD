% Example 12.4
clc
clear all      
% Starting estimates
    theta = [1.57; 0.5; 1.2];
% Compute the roots of the function     
    [coords, check] = NRn(@Ex_12_4, theta, 0.0001, 10);
    if check == 0
        'Convergence failed in NRn (Newton-Raphson)'
    else
        coordinates = coords
    end
