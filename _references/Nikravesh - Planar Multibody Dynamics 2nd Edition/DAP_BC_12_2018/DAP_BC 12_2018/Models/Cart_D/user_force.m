function user_force
    include_global

% Motor
    omega_max = 4*pi; T_max = 20;
    omega = abs(Bodies(2).p_d);
    T_motor = T_max*(1 - omega/omega_max);
    if T_motor > T_max
        T_motor = T_max;
    end
    Bodies(2).n = Bodies(2).n - T_motor;
    Bodies(1).n = Bodies(1).n + T_motor;

% Aerodynamic resistive force
    damp_aero = 10;
    x_d = Bodies(1).r_d(1);
    f_aero = damp_aero*x_d^2;
    Bodies(1).f(1) = Bodies(1).f(1) - f_aero;
    
end