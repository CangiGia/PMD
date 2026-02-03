function ff = Friction_B(mu_s, mu_d, mu_v, v_t, fnt, v, fn)
% Friction force based on Brown-McPhee model

% Viscous friction is not included
    vvt = v/v_t;
    ff = fn*(mu_d*tanh(4*vvt) + (mu_s - mu_d)*vvt/(0.25*vvt^2 + 0.75)^2) + ...
             mu_v*v*tanh(4*fn/fnt);
    
end
