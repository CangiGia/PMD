% disc_C
% Constraints for a disc joint  between a body and the flat ground

    Bi = Joints(Ji).iBindex;
    f = [(Bodies(Bi).r(2) - Joints(Ji).R)
         ((Bodies(Bi).r(1) - Joints(Ji).x0) + ...
           Joints(Ji).R*(Bodies(Bi).p - Joints(Ji).p0))];
