% revolute_C
% Evaluate the constraints for a revolute joint 

    Pi = Joints(Ji).iPindex;    Pj = Joints(Ji).jPindex;

    f = Points(Pi).rP - Points(Pj).rP;
    if Joints(Ji).fix == 1
        Bi = Joints(Ji).iBindex; Bj = Joints(Ji).jBindex;
        if Bi == 0
        f = [f
            (- Bodies(Bj).p - Joints(Ji).p0)];
        elseif Bj == 0
            f = [f
            (Bodies(Bi).p - Joints(Ji).p0)];
        else
            f = [f
            (Bodies(Bi).p - Bodies(Bj).p - Joints(Ji).p0)];
        end
    end