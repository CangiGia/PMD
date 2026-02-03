% translational_A
% r-h-s of acc. constraints for a translational joint

    Bi = Joints(Ji).iBindex;      Bj = Joints(Ji).jBindex;
    Pi = Joints(Ji).iPindex;      Pj = Joints(Ji).jPindex;
    ujd = Uvectors(Joints(Ji).jUindex).u_d;
    ujd_r = s_rot(ujd);        
        
    if Bi == 0
        f2 = 0;
    elseif Bj == 0
        f2 = 0;
    else
        f2 = ujd'*(Bodies(Bi).r - Bodies(Bj).r)*Bodies(Bi).p_d - ...
             2*ujd_r'*(Bodies(Bi).r_d - Bodies(Bj).r_d);
    end
        f  = [f2; 0];
 
    if Joints(Ji).fix == 1
        d    = Points(Pi).rP - Points(Pj).rP;
        d_d  = Points(Pi).rP_d - Points(Pj).rP_d;
        L = Joints(Ji).p0; u = d/L; u_d = d_d/L;
        f3 = - u_d'*d_d;
        if Bi == 0
            f3 = f3 + u'*s_rot(Points(Pj).sP_d)*Bodies(Bj).p_d;
        elseif Bj == 0
            f3 = f3 - u'*s_rot(Points(Pi).sP_d)*Bodies(Bi).p_d;
        else
            f3 = f3 - u'*(s_rot(Points(Pi).sP_d*Bodies(Bi).p_d - ...
                                Points(Pj).sP_d*Bodies(Bj).p_d));
        end
        f = [f; f3];
    end
    