function u_d = analysis(t, u)
% Solve the constrained equations of motion at time t with the standard
%   Lagrange multiplier method

    include_global
    
    u_to_Bodies; % unpack u into coordinate and velocity sub-arrays
    Update_Position; Update_Velocity;  
    h_a = Force_array(t); % array of applied forces
    if nConst == 0
        c_dd = M_inv_array.*h_a; % solve for accelerations
    else
    % Phi = Constraints(t);
        D = Jacobian; % Jacobian matrix
        rhsA = RHSAcc(t); % r-h-s of acc constraints (gamma)
        DMD = [diag(M_array) -D'
               D   zeros(nConst)];
        rhs = [h_a
               rhsA];
        sol = DMD\rhs;
        c_dd = sol(1:nB3);
        Lambda = sol(nB3+1:end);
    end
    for Bi=1:nB
        ir = Bodies(Bi).irc; i2 = ir + 1; i3 = i2 + 1;
        Bodies(Bi).r_dd = c_dd(ir:i2);
        Bodies(Bi).p_dd = c_dd(i3);
    end
    Bodies_to_u_d; % pack velocities and accelerations into ud
    
    num = num + 1; % number of function evaluations    
    if showtime == 1
    % Inform the user of the progress: 
    %   show the time once every 100 function evaluations
        if mod(t10, 100) == 0
            disp(t)
        end
            t10 = t10 + 1;
    end
