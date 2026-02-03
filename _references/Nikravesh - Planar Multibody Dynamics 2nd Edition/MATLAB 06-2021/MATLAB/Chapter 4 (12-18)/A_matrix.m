function A = A_matrix(phi)
% Computes rotational transformation matrix A
    cp = cos(phi); sp = sin(phi);
    A = [cp -sp; sp cp];
end
    