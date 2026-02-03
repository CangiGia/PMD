function ic_correct
% This function corrects initial conditions on the joint coordinates and
% velocities
    include_global
    global theta theta_d

% Coordinate correction
    flag = 0;
for n = 1:20
    ChainOpen(1, theta, theta_d, 0); Update_Position; % Update position entities
    CutJoint(1, theta, theta_d, 0);  % Evaluate Jacobian and constraints
    ff = sqrt(Phi'*Phi);  % Are the constraints violated?
    if ff < 1.0e-10
        flag = 1; break
    end
    delta_theta = -C'*((C*C')\Phi);  % Solve for corrections
    theta = theta + delta_theta;
end
    if flag == 0
        error(' Convergence failed in Newton-Raphson ');
    end
    
% Velocity correction
    CutJoint(2, theta, theta_d, 0);
    delta_thetad = -C'*((C*C')\(C*theta_d - rhsV)); % Solve for corrections
    theta_d = theta_d + delta_thetad;

% Report corrected results
    display(' ')
    display('Corrected joint coordinates')
    display(num2str(theta))
    display('Corrected joint velocities')
    display(num2str(theta_d))
    display(' ')
