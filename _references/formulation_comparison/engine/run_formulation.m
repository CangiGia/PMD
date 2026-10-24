function run_formulation(cfg)
%RUN_FORMULATION Run one DAP driver on one model and save its results.
%   The current folder must be the solver folder (containing dap.m). The
%   driver calls clear all, so the configuration is passed via appdata and
%   read again after the run. Output files: <formulation>_results.mat
%   and stats_<formulation>.txt in cfg.out_dir.
%
%   Both formulations save the same quantities so they can be compared:
%     q     - (x, y, phi) of every body, three columns per body
%     reac  - constraint force on every body (Fx, Fy, Mz), three columns
%             per body. BCF: D'*Lambda. JCF: M*c_dd - h_a.
%     lambda- raw Lagrange multipliers (BCF; JCF only with cut joints)
%     t, q_labels, q_kinds, reac_labels, reac_kinds
%   JCF integrates joint coordinates, so its q is rebuilt from the body
%   positions recorded by post.

    setappdata(0, 'cmp_cfg', cfg);
    setappdata(0, 'cmp_answers', struct('folder', cfg.folder, ...
        'ic', 'y', 'tf', cfg.tf, 'dt', cfg.dt));

    dap
    cfg = getappdata(0, 'cmp_cfg');  % dap runs clear all, which drops cfg

    % Replay the trajectory to fill the Lagrange multipliers (Lam).
    post

    if strcmp(cfg.formulation, 'BCF')
        q = uT(:, 1:nB3);
        reac = bcf_reactions(Jac, Lam, nConst, nB3);
    else
        q = jcf_body_coords(r, p);
        reac = jcf_reactions(uT, T, theta_dd);
    end
    lambda = Lam;
    t = T;

    % dap.m evaluates Phi after ChainOpen only, without Update_Position,
    % which leaves some points stale. Recompute the residual with the
    % positions updated, as ic_correct does.
    if strcmp(cfg.formulation, 'JCF') && nConst ~= 0
        phi_norm = jcf_phi_updated(uT, T);
    end

    phi_max = max(phi_norm);
    phi_end = phi_norm(end);
    lam_max = 0;
    if ~isempty(Lam)
        lam_max = max(abs(Lam(:)));
    end

    [q_labels, q_kinds] = coord_labels(size(q, 2));
    [reac_labels, reac_kinds] = reaction_labels(size(reac, 2));

    save(fullfile(cfg.out_dir, [cfg.formulation '_results.mat']), ...
        'q', 'reac', 'lambda', 't', 'q_labels', 'q_kinds', ...
        'reac_labels', 'reac_kinds');

    if nConst ~= 0
        rank_txt = sprintf('%d of %d', redund, nConst);
    else
        rank_txt = 'not applicable (no constraints)';
    end

    fid = fopen(fullfile(cfg.out_dir, ['stats_' cfg.formulation '.txt']), 'w');
    fprintf(fid, 'Model                       : %s\n', cfg.model);
    fprintf(fid, 'Formulation                 : %s\n', cfg.formulation);
    fprintf(fid, 'Model folder                : %s\n', cfg.folder);
    fprintf(fid, 'Solver folder               : %s\n', cfg.solver_dir);
    fprintf(fid, 'Final time [s]              : %g\n', cfg.tf);
    fprintf(fid, 'Reporting time-step [s]     : %g\n', cfg.dt);
    fprintf(fid, 'Time samples                : %d\n', numel(t));
    fprintf(fid, 'Generalized coordinates     : %d\n', size(q, 2));
    fprintf(fid, 'Number of constraints       : %d\n', nConst);
    fprintf(fid, 'Constraint rank             : %s\n', rank_txt);
    fprintf(fid, 'Function evaluations        : %d\n', num);
    fprintf(fid, 'Simulation time [s]         : %.6g\n', sim_time);
    fprintf(fid, 'Max constraint violation    : %.6g\n', phi_max);
    fprintf(fid, 'Final constraint violation  : %.6g\n', phi_end);
    fprintf(fid, 'Max |lambda|                : %.6g\n', lam_max);
    fprintf(fid, 'Generated                   : %s\n', datestr(now));
    fclose(fid);

    fprintf('Saved %s results for %s in %s\n', cfg.formulation, ...
        cfg.model, cfg.out_dir);
end

function phi = jcf_phi_updated(uT, T)
%JCF_PHI_UPDATED Position-constraint norm ||Phi|| at each reporting time.
%   Positions are updated (Update_Position) after ChainOpen so that every
%   point used by CutJoint is current. Requires the globals set by dap.

    include_global
    phi = zeros(numel(T), 1);
    for k = 1:numel(T)
        th = uT(k, 1:nJC)';
        th_d = uT(k, nJC+1:nJC2)';
        ChainOpen(1, th, th_d, T(k));
        Update_Position;
        CutJoint(1, th, th_d, T(k));
        phi(k) = norm(Phi);
    end
end

function reac = bcf_reactions(Jac, Lam, nConst, nB3)
%BCF_REACTIONS Constraint force on each body, D'*Lambda, at each time.
%   Jac is nt x nConst x nB3 and Lam is nt x nConst, as saved by post.

    nt = size(Jac, 1);
    reac = zeros(nt, nB3);
    if nConst == 0
        return
    end
    for k = 1:nt
        D = reshape(Jac(k, :, :), nConst, nB3);
        reac(k, :) = (D' * Lam(k, :)')';
    end
end

function q = jcf_body_coords(r, p)
%JCF_BODY_COORDS Body (x, y, phi) per body from the JCF post output.
%   r is nt x nB x 2 and p is nt x nB. The columns follow the BCF order:
%   body 1 x, y, phi, body 2 x, y, phi, and so on.

    nt = size(p, 1);
    nb = size(p, 2);
    q = zeros(nt, 3 * nb);
    for b = 1:nb
        q(:, 3*b-2:3*b-1) = reshape(r(:, b, :), nt, 2);
        q(:, 3*b) = p(:, b);
    end
end

function reac = jcf_reactions(uT, T, theta_dd)
%JCF_REACTIONS Constraint force on each body, M*c_dd - h_a, at each time.
%   Follows analysis.m: ChainOpen(2), Update_Velocity, ChainOpen(3),
%   Force_array. c_dd = B*theta_dd + Bd. The open chain has no Lagrange
%   multipliers, so this is the only source of the constraint forces.
%   Requires the globals set by dap and theta_dd from post.

    include_global
    nt = numel(T);
    reac = zeros(nt, nB3);
    for k = 1:nt
        th = uT(k, 1:nJC)';
        th_d = uT(k, nJC+1:nJC2)';
        ChainOpen(2, th, th_d, T(k));
        Update_Position;
        B = Bmat;
        c_d = B * th_d;
        Update_Velocity;
        ChainOpen(3, th, th_d, T(k));
        h_a = Force_array(T(k));
        c_dd = B * theta_dd(k, :)' + Bd;
        reac(k, :) = (M_array .* c_dd - h_a)';
    end
end
