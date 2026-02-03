% Filename: Exaple 13_9: Two pendulums impacting each other
    clc; clear all

    global M_diag M h_a L R R2 e

% Data
    L = 0.1; R = 0.01; R2 = 2*R; e = 1; g = 9.81;
    m1 = 0.1; m2 = m1; J1 = 0.5*m1*R^2; J2 = J1;
    M_diag = [m1 m1 J1 m2 m2 J2]; M = diag(M_diag);
    h_a = [0 -m1*g 0 0 -m2*g 0]';
    
% Initial coordinates
    p1 = pi/3; r1 = [R2; 0] + L*[sin(p1); -cos(p1)];
    p2 = 0; r2 = [0; -L];
% Integration array
    u = [r1' p1 r2' p2  zeros(1,6)]';
% Time span
    Tstart = 0; Tfinal = 2; delT = 0.01;
% Integration options
    options = odeset('Event', @Ex_13_9_2, 'RelTol', 1e-6, 'AbsTol', 1e-9);
    
    TT = []; uTT = []; % storage arrays for segmented results
    
while Tstart < Tfinal
    Tspan = [Tstart Tfinal]; % report results at every time step
%     Tspan = [Tstart: delT: Tfinal];
    [T, uT, Te, uTe, ie] = ode45(@Ex_13_9_1, Tspan, u, options);
    TT = [TT; T]; uTT = [uTT; uT];
    if T(end)+delT < Tfinal
        vp = Ex_13_9_3(uTe(end,:)', M, R, L, e);
        u = [uTe(end,1:6)';vp];
        Tstart = Te(end);
    else
        break
    end
end
%     subplot(4,1,1)
%     plot(TT, uTT(:,3),'k')
%     subplot(4,1,2)
%     plot(TT, uTT(:,6),'k')
%     subplot(4,1,3)
%     plot(TT, uTT(:,9),'k')
%     subplot(4,1,4)
%     plot(TT, uTT(:,12),'k')
    
    Ex_13_9_4
    
%     nT = length(TT);
%     dTs = []; TTT = [];
%     dTs = TT(2:end) - TT(1:end-1);
% 
%     figure
%     subplot(2,1,1)
%     plot(TT, uTT(:,9),'k')
%     subplot(2,1,2)
%     plot(TT(1:end-1),dTs)
    