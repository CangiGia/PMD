function [coords, check] = NRn(fun, q, tol, iter_max)
% Newton-Raphson iterative process
    
    check = 0;
for i = 1:iter_max 
        [Phi, D] = fun(q);         % constraints and Jacobian
        ff = norm(Phi);            % are the constraints violated?
    if ff < tol
        check = 1;                 % a solution found!!!
        break
    end
        delta_q = -D\Phi;          % solve for corrections
        q = q + delta_q;           % correct the estimates
end
    coords = q;
