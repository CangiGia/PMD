function ff = Friction_B(mu_s, mu_d, mu_v, v_t, fnt, v, fN)
% Friction force based on Brown-McPhee model

% Viscous friction is included
    vr = v/v_t;
    ff = fN*(mu_d*tanh(4*vr) + (mu_s - mu_d)*vr/(0.25*vr^2 + 0.75)^2) + ...
             mu_v*v*tanh(4*fN/fnt);
end
