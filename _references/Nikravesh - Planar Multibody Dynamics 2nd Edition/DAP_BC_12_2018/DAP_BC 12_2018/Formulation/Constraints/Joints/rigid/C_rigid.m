% rigid_C
% Constraints for a Rigid joint

    Bi = Joints(Ji).iBindex;      Bj = Joints(Ji).jBindex;
    
    if Bi == 0
        f = [ -(Bodies(Bj).r + Bodies(Bj).A*Joints(Ji).d0)
              -Bodies(Bj).p - Joints(Ji).p0];
    elseif Bj == 0
        f = [Bodies(Bi).r - Joints(Ji).d0
             Bodies(Bi).p - Joints(Ji).p0];
    else
        f = [Bodies(Bi).r - (Bodies(Bj).r + Bodies(Bj).A*Joints(Ji).d0)
             Bodies(Bi).p - Bodies(Bj).p - Joints(Ji).p0];
    end

