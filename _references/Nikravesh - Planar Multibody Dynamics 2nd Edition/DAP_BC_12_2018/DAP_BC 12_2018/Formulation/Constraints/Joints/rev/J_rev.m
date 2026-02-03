% revolute_J
% Jacobian sub-matrices for a revolute joint

Pi = Joints(Ji).iPindex;    Pj = Joints(Ji).jPindex;

    Di = [ eye(2)  Points(Pi).sP_r];
    Dj = [-eye(2) -Points(Pj).sP_r];
    
    if Joints(Ji).fix == 1
        Di = [Di
              0  0  1];
        Dj = [Dj
              0  0 -1];
    end