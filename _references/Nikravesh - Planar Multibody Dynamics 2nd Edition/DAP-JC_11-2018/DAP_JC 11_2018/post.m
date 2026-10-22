% Save joint  coordinates, velocities, and accelerations
%     and the coordinates and velocity of all points
%     Jacobian matrix and Lagrange multipliers 
%     at every reporting time interval

    nt  = size(T,1);          % number of time steps
    theta    = zeros(nt,nJC); % joint coordinates
    theta_d  = zeros(nt,nJC); % joint velocities
    theta_dd = zeros(nt,nJC); % joint accelerations
    r   = zeros(nt,nB,2);     % translational coordinates
    rd  = zeros(nt,nB,2);     % translational velocities
%     rdd = zeros(nt,nB,2);   % translational acceleration
    p   = zeros(nt,nB);       % rotational coordinate
    pd  = zeros(nt,nB);       % angular velocity
%     pdd = zeros(nt,nB);     % angular acceleration
    rP  = zeros(nt,nP,2);     % coordinates of points
    rPd = zeros(nt,nP,2);     % velocity of points
    Jac = zeros(nt,nConst,nJC); % Jacobian matrix
    Lam = zeros(nt,nConst);   % Lagrange multipliers

    showtime = 0;
    for i=1:nt
        t = T(i); u = uT(i,:)'; % u_to_Bodies;
        theta(i,:) = u(1:nJC);
        theta_d(i,:) = u(nJC+1:2*nJC);        
        u_d = analysis(t, u);
        theta_dd(i,:) = u_d(nJC+1:2*nJC)';
        for Bi = 1:nB;
            r(i,Bi,:)   = Bodies(Bi).r;
            p(i,Bi)     = Bodies(Bi).p;
            rd(i,Bi,:)  = Bodies(Bi).r_d;
            pd(i,Bi)    = Bodies(Bi).p_d;
%             rdd(i,Bi,:) = Bodies(Bi).r_dd;
%             pdd(i,Bi)   = Bodies(Bi).p_dd;
        end
        for j=1:nP
            rP(i,j,:)   = Points(j).rP;
            rPd(i,j,:)  = Points(j).rP_d;
        end
        if nConst > 0
            Jac(i,:,:) = C;
            Lam(i,:) = Lambda';
        end
    end
