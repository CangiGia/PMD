% Examplle 12.6 (a continuation of Example 12.1)
% Execute the program from Example 12.1
    Example_12_1;
% Additional constant data
    m1 = 1.0; m2 = 1.5; m3 = 10.0;
    J1 = 0.001; J2 = 0.002; J3 = 1.0;
    f3 = 10;
% Initialize arrays
    Lag = zeros(nsteps,8); Torque1 = zeros(nsteps,1); Torque2 = zeros(nsteps,1);
for i = 1:nsteps
% Retrieve kinematic results
    theta    = thetas(i,:);
    theta_d  = thetas_d(i,:);
    theta_dd = thetas_dd(i,:);
% Unit vectors
    c1 = cos(theta(1)); s1 = sin(theta(1)); 
    c2 = cos(theta(2)); s2 = sin(theta(2));
    u1  = [ c1; s1]; u2  = [ c2; s2]; u0  = [ 1; 0];
    u1r = [-s1; c1]; u2r = [-s2; c2]; u0r = [ 0; 1];
% Velocity of mass centers
    rG1_d = 0.5*L1*u1r*theta_d(1);
    rG2_d = 2*rG1_d + 0.5*L2*u2r*theta_d(2);
    rG3_d = u0*theta_d(3);
% Acceleration of mass centers
    rG1_dd = 0.5*L1*(u1r*theta_dd(1) - u1*theta_d(1)^2);
    rG2_dd = 2*rG1_dd + 0.5*L2*(u2r*theta_dd(2) - u2*theta_d(2)^2);
    rG3_dd = u0*theta_dd(3);
% Right-hand side
    rhs = [m1*rG1_dd
           J1*theta_dd(1) 
           m2*rG2_dd
           J2*theta_dd(2)
           m3*rG3_dd - f3*u0
           0];
%% First method
% Moment arms (rotated)
    sO1r = -0.5*L1*u1r; sA1r = -sO1r;
    sA2r = -0.5*L2*u2r; sB2r = -sA2r;
% Coefficient matrix (Jacobian transpose)
    I22 = eye(2); Z22 = zeros(2); Z12 = zeros(1,2);
    Dt = [I22    I22    Z22    Z22       Z12'
          sO1r'  sA1r'  Z12    Z12       1
          Z22   -I22    I22    Z22       Z12'
          Z12   -sA2r'  sB2r'  Z12       0
          Z22    Z22   -I22    Z12' u0r  Z12'
          Z12    Z12    Z12    1 0       0  ];

% Solve equations for Lagrange multipliers
    sol = Dt\rhs;
    Lag(i,:) = sol(1:8)';
    Torque1(i) = sol(9);
%% Second method
% Construct array of velocities
    c_d = [rG1_d
           theta_d(1)
           rG2_d
           theta_d(2)
           rG3_d
           0 ];
% Solve equation for one Lagrange multiplier
    Torque2(i) = c_d'*rhs/theta_d(1);
end
subplot(1,2,1)
plot(T, Torque1)
subplot(1,2,2)
plot(T, Torque2)

