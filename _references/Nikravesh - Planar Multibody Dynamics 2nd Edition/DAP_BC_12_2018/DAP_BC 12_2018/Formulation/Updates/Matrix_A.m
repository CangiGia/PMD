function A = Matrix_A(p)
% This function computes the rotational transformation matrix A
    c =cos(p); s = sin(p);
        A = [c -s
             s  c];
     