% SDA_ptp
% Point-to-point spring-damper-actuator
    Pi     = Forces(Fi).iPindex;    Pj     = Forces(Fi).jPindex;
    Bi     = Forces(Fi).iBindex ;   Bj     = Forces(Fi).jBindex;
    d      = Points(Pi).rP - Points(Pj).rP;        
    d_dot  = Points(Pi).rP_d - Points(Pj).rP_d;   
    L     = sqrt(d'*d);
    L_dot = d'*d_dot/L;
    del    = L - Forces(Fi).L0;
    u      = d/L;
    
    f = Forces(Fi).k*del + Forces(Fi).dc*L_dot + Forces(Fi).f_a;
    fi = f*u;
    if Bi ~= 0
        Bodies(Bi).f = Bodies(Bi).f - fi;
        Bodies(Bi).n = Bodies(Bi).n - Points(Pi).sP_r'*fi;
    end
    if Bj ~= 0
        Bodies(Bj).f = Bodies(Bj).f + fi;
        Bodies(Bj).n = Bodies(Bj).n + Points(Pj).sP_r'*fi;
    end
