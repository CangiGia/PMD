function A = A_matrix(phi)
% This function computes rotational matrix A
    cp = cos(phi); sp = sin(phi);
    A = [cp -sp; sp cp];
end
    