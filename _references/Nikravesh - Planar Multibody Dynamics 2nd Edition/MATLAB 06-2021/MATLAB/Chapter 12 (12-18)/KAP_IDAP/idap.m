%% idap.m
% Kinematic and Inverse Dynamic Analysis of a four-bar mechanism 
%   (Appended Constraint method)
% One coupler point and one follower point 
% The user must supply the following data:
%   Constant lengths: L = [L1 L2 L3]
%       RQ = [L0x L0y] (ground), L1(crank), L2 (coupler), L3 (follower)
%   Secondary points on links 1, 2 and 3: 
%       LP = [LP1 LP2 LP3] (lengths)
%       beta = [beta1 beta2 beta3] (angles in degrees)
%   Angles and angular velocity
%       [theta1 (exact) theta2 and theta3 (estimates)] (in degrees)
%   Angular velocity of the crank
%       omega1 (in rad/sec)
%   Parameters for plot (animation) window
%   Masses m = [m1 m2 m3]; Moments of inertia J = [J1 J2 J3]
%   Position of mass centers LG = [LG1 :G2 LG3]
%   Angle of mass centers gamma = [gamma1 gamma2 gamma3]
%   Gravitational constant g
%   Directional unit vector for gravity ug = [ux; uy]
%%
    clc; clear all
    global RQ R_A R_B R_P1 R_P2 R_P3 P1_yes P2_yes P3_yes axis_limits
    global L theta1
    addpath Formulation Models
% User supplied data 
    folder = input(' Which folder contains the model? ', 's');
    addpath(fullfile(pwd,'Models',folder))
    data 
    RQ = RQ(:); theta = theta(:); beta = beta(:);
% Transform angles to radians
    theta = theta*pi/180; beta = beta*pi/180; theta1 = theta(1);
    del_theta1 = 2*sign(omega1);   % Increments of 2 degrees for the crank
    del_theta1 = del_theta1*pi/180; theta1_0 = theta(1);
% Grashof condition
    L0 = norm(RQ); Lengths = sort([L0; L']);
    sl = Lengths(1) + Lengths(4); pq = Lengths(2) + Lengths(3);
    if sl <= pq         disp(' The four-bar is Grashof')
    else                disp(' The four-bar is non-Grashof')
    end
% Check for points P and C
    P1_yes = 1; if LP(1) == 0; P1_yes = 0; end
    P2_yes = 1; if LP(2) == 0; P2_yes = 0; end
    P3_yes = 1; if LP(3) == 0; P3_yes = 0; end
% Initialize
    nsteps = (2*pi/abs(del_theta1) + 1); del_time = (2*pi/abs(omega1))/nsteps;
    T = zeros(nsteps,1);     time = 0;
    thetas = zeros(nsteps,3); omegas = zeros(nsteps,3); alphas = zeros(nsteps,3); 
    
    RA = zeros(nsteps,2); RB = zeros(nsteps,2); 
    VA = zeros(nsteps,2); VB = zeros(nsteps,2); 
    AA = zeros(nsteps,2); AB = zeros(nsteps,2); 
    
    RP1 = zeros(nsteps,2); RP2 = zeros(nsteps,2); RP3 = zeros(nsteps,2); 
    VP1 = zeros(nsteps,2); VP2 = zeros(nsteps,2); VP3 = zeros(nsteps,2); 
    AP1 = zeros(nsteps,2); AP2 = zeros(nsteps,2); AP3 = zeros(nsteps,2); 

    RF = zeros(nsteps,8); Torque1 = zeros(nsteps,1); Torque2 = zeros(nsteps,1);

    RG1 = zeros(nsteps,2); RG2 = zeros(nsteps,2); RG3 = zeros(nsteps,2); 
    VG1 = zeros(nsteps,2); VG2 = zeros(nsteps,2); VG3 = zeros(nsteps,2); 
    AG1 = zeros(nsteps,2); AG2 = zeros(nsteps,2); AG3 = zeros(nsteps,2); 
% Mass and array of weights
    mass = [m(1) m(1) J(1) m(2) m(2) J(2) m(3) m(3) J(3)]';
    weight = [m(1)*g*ug; 0; m(2)*g*ug; 0; m(3)*g*ug; 0];
for i = 1:nsteps
%% Start of Kinematic Analysis% Position analysis
%
    [theta, check] = NRn(@constraints, theta, 0.0001, 10);
    if check == 0
        del_theta1 = -del_theta1; % Change the direction of the crank's rotation
        theta1 = theta1 + del_theta1;
        [theta, check] = NRn(@constraints, theta, 0.0001, 10);
        if check == 0
            'Convergence failed in NRn (Newton-Raphson)'
            break
        end
    end
% Velocity analysis: Solve for unknown velocities
    [u, ur] = uVectors(theta);  % update unit vectors
    D  = [L(1)*ur{1} L(2)*ur{2} -L(3)*ur{3}
          1            0             0        ];  %  update Jacobian
    rhsv = [0; 0; omega1]; % Right-hand-side of velocity constraints
    omega = D\rhsv;
% Acceleration analysis: Solve for unknown accelerations
    rhsa = [L(1)*u{1}*omega(1)^2 + L(2)*u{2}*omega(2)^2 - L(3)*u{3}*omega(3)^2
            0];            % Right-hand-side of acceleration constraints
    alpha = D\rhsa;  
% Points A, B, P1, P2, and P3
    Pees
% Save results
    T(i) = time;
    thetas(i,:) = theta; omegas(i,:) = omega'; alphas(i,:) = alpha';
    
%% --------- End of kinematic analysis -----------------  
% 
%% Start of Inverse Dynamic Analysis
%
% Transform angles to radians
    gamma = gamma(:); gamma = gamma*pi/180;
% Coordinates, velocity and acceleration of mass centers 
    Gees
% Arrays of body velocities and accelerations
    c_d  = [VG1(i,:) omegas(i,1) VG2(i,:) omegas(i,2) VG3(i,:) omegas(i,3)]';
    c_dd = [AG1(i,:) alphas(i,1) AG2(i,:) alphas(i,2) AG3(i,:) alphas(i,3)]';
% Array of applied forces
    h_a = weight;
% Evaluate moment arms and the 9x9 Jacobian
    moms
%  Include user-supplied forces to the array h_a
    user_force
% Solve for the torque of the motor
        Torque1(i) = c_d'*(mass.*c_dd - h_a)/omega1;
% Solve for the torque of the motor and reaction forces
        Lag = D'\(mass.*c_dd - h_a);
        RF(i,:) = Lag(1:8)';
        Torque2(i) = Lag(9);
%% --------- End of inverse dynamic analysis -----------------  
%
% Increment the crank angle and time
    theta1 = theta1 + del_theta1;
    time = time + del_time;
end
    hold off
%     subplot(1,2,1)
    plot(T,Torque1,'k')
%     subplot(1,2,2)
%     plot(T,Torque2,'r')
