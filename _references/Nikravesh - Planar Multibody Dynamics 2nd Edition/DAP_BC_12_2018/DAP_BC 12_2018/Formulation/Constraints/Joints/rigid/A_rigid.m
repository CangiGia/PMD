% rigid_A
% r-h-s of acc. constraints for a Rigid joint

    Bj = Joints(Ji).jBindex;
    
    f = [0; 0; 0];
    if Bj ~= 0
        f = [-Bodies(Bj).A*Joints(Ji).d0*Bodies(Bj).p_d^2; 0];

    end