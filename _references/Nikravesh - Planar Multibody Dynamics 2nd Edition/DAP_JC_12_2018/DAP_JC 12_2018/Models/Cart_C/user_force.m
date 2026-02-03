function user_force
    include_global

% Motor
    omega_max = 4*pi; T_max = 20;
    omega = abs(Bodies(1).p_d);
    T_motor = T_max*(1 - omega/omega_max);
    if T_motor > T_max
        T_motor = T_max;
    end
    Bodies(1).n = Bodies(1).n - T_motor;
    Bodies(2).n = Bodies(2).n + T_motor;    
end