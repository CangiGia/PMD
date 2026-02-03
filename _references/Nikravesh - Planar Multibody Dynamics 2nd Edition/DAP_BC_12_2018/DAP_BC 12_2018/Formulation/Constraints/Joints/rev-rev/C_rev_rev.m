% revolute_revolute_C
% Evaluates the constraint for a revolute_revolute joint 

    Pi = Joints(Ji).iPindex;  Pj = Joints(Ji).jPindex;
    d  = Points(Pi).rP - Points(Pj).rP;
    L = Joints(Ji).L; u = d/L;
        f = (u'*d - L)/2;
    