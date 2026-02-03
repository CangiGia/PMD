function u_d = analysis(t, u)
% Solve the constrained equations of motion at time t with the standard
%   Lagrange multiplier method

    include_global
    theta = u(1:nJC); theta_d = u(nJC+1:nJC2);
    ChainOpen(2, theta, theta_d, 0);

    B = Bmat; 
    c_d = B*theta_d; 
    Update_Velocity; 
    ChainOpen(3, theta, theta_d, 0);
    
    h_a = Force_array(t);  % array of applied forces
    h_joint = B'*(h_a - M_array.*Bd);
    M_joint = B'*diag(M_array)*B;
    if nConst == 0
        theta_dd = M_joint\h_joint; % solve for accelerations
    else
    % CutJoint(1, theta, theta_d, t); % Jacobian and constraints
        CutJoint(3, theta, theta_d, t); % Jacobian and r-h-s of acc.
        CMC = [M_joint -C'
               C   zeros(nConst)];
        rhs = [h_joint
               rhsA];
        sol = CMC\rhs;
        theta_dd = sol(1:nJC);
        Lambda = sol(nJC+1:end);
    end

    u_d = [theta_d; theta_dd];
    
    num = num + 1; % number of function evaluations    
    if showtime == 1
    % Inform the user of the progress: 
    %   show the time once every 100 function evaluations
        if mod(t10, 100) == 0
            disp(t)
        end
            t10 = t10 + 1;
    end
