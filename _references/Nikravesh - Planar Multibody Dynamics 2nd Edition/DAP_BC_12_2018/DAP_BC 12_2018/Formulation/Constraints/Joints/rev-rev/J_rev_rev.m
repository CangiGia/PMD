% rev_rev_J
% Jacobian sub-matrices for a revolute_revolute joint

    Pi = Joints(Ji).iPindex;  Pj = Joints(Ji).jPindex;
    d  = Points(Pi).rP - Points(Pj).rP;
    L = Joints(Ji).L; u = d/L;
        Di = [ u'  u'*Points(Pi).sP_r];
        Dj = [-u' -u'*Points(Pj).sP_r];
