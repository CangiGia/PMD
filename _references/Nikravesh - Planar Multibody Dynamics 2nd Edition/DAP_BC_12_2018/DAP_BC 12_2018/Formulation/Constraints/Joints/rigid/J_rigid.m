% rigid_J
% Jacobian sub-matrix for a Rigid joint  

    Bj = Joints(Ji).jBindex;
    
    Di = eye(3);
    if Bj ~= 0
        Dj = [-eye(2) -s_rot(Bodies(Bj).A*Joints(Ji).d0)
               0  0   -1];
    end
