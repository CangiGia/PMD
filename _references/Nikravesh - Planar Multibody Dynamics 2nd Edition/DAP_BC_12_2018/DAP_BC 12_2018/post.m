% Post Processor
% Compute body accelerations, Lagrange mulyipliers, coordinates and 
%    velocity of all points, kinetic and potential energies,
%    at every reporting time interval

    nt  = size(T,1);        % number of time steps
    r   = zeros(nt,nB,2);   % translational coordinates
    rd  = zeros(nt,nB,2);   % translational velocities
    rdd = zeros(nt,nB,2);   % translational acceleration
    p   = zeros(nt,nB);     % rotational coordinate
    pd  = zeros(nt,nB);     % angular velocity
    pdd = zeros(nt,nB);     % angular acceleration
    rP  = zeros(nt,nP,2);   % coordinates of points
    rPd = zeros(nt,nP,2);   % velocity of points
    Jac = zeros(nt,nConst,nB3); % Jacobian matrix
    Lam = zeros(nt,nConst); % Lagrange multipliers
    eng = zeros(nt,3);      % Energy (kinetic, potential, total)

    showtime = 0;
for i=1:nt
    t = T(i); u = uT(i,:)'; u_to_Bodies;
    analysis(t, u);
    for Bi = 1:nB;
        r(i,Bi,:)   = Bodies(Bi).r;    p(i,Bi)   = Bodies(Bi).p;
        rd(i,Bi,:)  = Bodies(Bi).r_d;  pd(i,Bi)  = Bodies(Bi).p_d;
        rdd(i,Bi,:) = Bodies(Bi).r_dd; pdd(i,Bi) = Bodies(Bi).p_dd;
    end
    for j=1:nP
        rP(i,j,:)   = Points(j).rP; rPd(i,j,:)  = Points(j).rP_d;
    end
    if nConst > 0
        Jac(i,:,:) = D; Lam(i,:) = Lambda';
    end
    % Compute kinetic and potential energies
    kin = uT(i,nB3+1:end)*(M_array.*uT(i,nB3+1:end)')/2;
    potential = 0;
    for Fi = 1:nF
        switch (Forces(Fi).type);
            case {'weight'}
                for Bi=1:nB
                    potential = potential - Bodies(Bi).wgt'*Bodies(Bi).r;
                end
            case {'ptp'}
                SDA_ptp
                potential = potential + 0.5*Forces(Fi).k*del^2;
            case {'rot-sda'}
                SDA_rot
            case {'flocal'}
                Bi = Forces(Fi).iBindex;
                Bodies(Bi).f = Bodies(Bi).f + Bodies(Bi).A*Forces(Fi).flocal; 
            case {'f'}
                Bi = Forces(Fi).iBindex;
                Bodies(Bi).f = Bodies(Bi).f + Forces(Fi).f; 
            case {'T'}
                Bi = Forces(Fi).iBindex;
                Bodies(Bi).n = Bodies(Bi).n + Forces(Fi).T; 
            case {'user'}
                user_force
        end
    end
    tot = kin + potential;
    eng(i,:) = [kin, potential, tot];
end
