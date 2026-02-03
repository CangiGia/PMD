function [u, ur] = uVectors(theta)
% compute unit and rotated unit vectors
    u = {}; ur = {};
    for i = 1:3
        c = cos(theta(i)); s = sin(theta(i));
        u{i}  = [c; s];
        ur{i} = [-s; c];
    end
