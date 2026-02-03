% SDA_rot
% Rotational spring-damper-actuator

    Bi = Forces(Fi).iBindex; Bj = Forces(Fi).jBindex;
    
    if Bi == 0
        theta   = -Bodies(Bj).p;
        theta_d = -Bodies(Bj).p_d; 
        T = Forces(Fi).k*(theta - Forces(Fi).theta0) + ...
            Forces(Fi).dc*theta_d + Forces(Fi).T_a;
        Bodies(Bj).n = Bodies(Bj).n + T;
    elseif Bj == 0
        theta   = Bodies(Bi).p;
        theta_d = Bodies(Bi).p_d; 
        T = Forces(Fi).k*(theta - Forces(Fi).theta0) + ...
            Forces(Fi).dc*theta_d + Forces(Fi).T_a;
        Bodies(Bi).n = Bodies(Bi).n - T;
    else
        theta   = Bodies(Bi).p - Bodies(Bj).p;
        theta_d = Bodies(Bi).p_d - Bodies(Bj).p_d; 
        T = Forces(Fi).k*(theta - Forces(Fi).theta0) + ...
            Forces(Fi).dc*theta_d + Forces(Fi).T_a;
        Bodies(Bi).n = Bodies(Bi).n - T;
        Bodies(Bj).n = Bodies(Bj).n + T;
    end