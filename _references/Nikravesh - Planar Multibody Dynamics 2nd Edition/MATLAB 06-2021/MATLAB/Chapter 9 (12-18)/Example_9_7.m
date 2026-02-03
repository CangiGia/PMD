% Example 9.7 (Initial Conditions: slider-crank mechanism)
    clc; clear all

    global r_O L1 L3
% Constant data
    r_O = [0.15; 0.33]; 
    L1 = 0.208; L3 = 0.195;
% Crank angle and angular velocity
    theta1 = 345*pi/180;
    theta1_d = pi;
% estimate for the other two coordinates
    theta23 = [0.2; 0.4];
% fsolve's parameters
    options = optimset('display', 'off');
    coords = fsolve(@s_c, theta23, options, theta1);
    angles_deg = coords*180/pi;
    theta2 = coords(1)
    theta3 = coords(2)
% Unit vectors
    u1 = [cos(theta1); sin(theta1)];
    u3 = [cos(theta2); sin(theta2)]; u3r = s_rot(u3);
% Jacibian matrix
    C = [L1*s_rot(u1) (L3*u3 - theta3*u3r) -u3
         1  0  0];
% Right-hand-side of velocity constraints
    rhsv = [0; 0; theta1_d];
% Solve for unknown velocities
    velocities = C\rhsv;
    theta2_d = velocities(2)
    theta3_d = velocities(3)
