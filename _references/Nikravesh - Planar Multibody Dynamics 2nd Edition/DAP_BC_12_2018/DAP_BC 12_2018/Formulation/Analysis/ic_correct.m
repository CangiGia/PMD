function ic_correct
% This function corrects initial conditions on the body coordinates and
% velocities
    include_global
    global Phi D rhs % this global is for debugging purpose

% Coordinate correction
    flag = 0;
for n = 1:20
    Update_Position;      % Update position entities
    Phi = Constraints(0); % Evaluate constraints
    D   = Jacobian;       % Evaluate Jacobian
    ff = sqrt(Phi'*Phi);  % Are the constraints violated?
    if ff < 1.0e-10
        flag = 1; break
    end
        delta_c = -D'*((D*D')\Phi);    % Solve for corrections
    for Bi = 1:nB         % Correct estimates
        ir = 1 + (Bi - 1)*3; 
        Bodies(Bi).r = Bodies(Bi).r + delta_c(ir:ir+1);
        Bodies(Bi).p = Bodies(Bi).p + delta_c(ir+2);
    end
end
    if flag == 0
        error(' Convergence failed in Newton-Raphson ');
    end

% Velocity correction
    for Bi = 1:nB     % Move velocities to an arbitrary array Phi
        ir = 1 + (Bi - 1)*3;
        Phi(ir:ir+2,1) = [Bodies(Bi).r_d; Bodies(Bi).p_d];
    end
        rhs = RHSVel(0);
        delta_v = -D'*((D*D')\(D*Phi - rhs)); % Compute corrections
    for Bi = 1:nB     % Move corrected velocities to sub-arrays
        ir  = 1 + (Bi - 1)*3;
        Bodies(Bi).r_d = Bodies(Bi).r_d + delta_v(ir:ir+1);
        Bodies(Bi).p_d = Bodies(Bi).p_d + delta_v(ir+2);
    end

% Report corrected coordinates and velocities
    coords = zeros(nB,3); vels = zeros(nB,3);
    for Bi = 1:nB
        coords(Bi,:) = [Bodies(Bi).r'  Bodies(Bi).p];
        vels(Bi,:) = [Bodies(Bi).r_d'  Bodies(Bi).p_d];
    end
    display(' ')
    display('Corrected coordinates')
    display(' x           y           phi')
    display(num2str(coords))
    display('Corrected velocities')
    display(' x-dot       y-dot       phi-dot')
    display(num2str(vels))
    display(' ')
