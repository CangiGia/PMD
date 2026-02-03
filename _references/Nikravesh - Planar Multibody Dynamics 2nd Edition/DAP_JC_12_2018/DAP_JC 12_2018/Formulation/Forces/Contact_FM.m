function fn = Contact_FM(del, deld, deld0, K, e)
% Contact force model Flores-Machado-Silva-Martins

    fn = K*del^1.5*(1 + 8*(1 - e)*deld/(5*e*deld0));
    
end
