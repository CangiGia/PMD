function ff = Friction_A(mu_s, mu_d, v_s, p, k_t, v, fN)
% Friction force based on Anderson et al. model

% Viscous friction is not included
    ff = fN*(mu_d + (mu_s - mu_d)*exp(-(abs(v)/v_s)^p))*tanh(k_t*v);
    
end
