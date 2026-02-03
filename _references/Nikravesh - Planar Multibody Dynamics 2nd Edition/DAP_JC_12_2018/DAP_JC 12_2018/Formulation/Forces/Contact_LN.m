function fn = Contact_LN(del, deld, deld0, K, e)
% Contact force model Lankarani-Nikravesh

    fn = K*del^1.5*(1 + 3*(1 - e^2)*deld/(4*deld0));
    
end
