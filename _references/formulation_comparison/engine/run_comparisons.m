function run_comparisons(only)
%RUN_COMPARISONS Run the BCF and JCF simulations for every paired model.
%   run_comparisons() runs all paired models in both formulations.
%   run_comparisons('MODEL') runs a single model (folder name from the
%   BCF column of the table below). Results are written to
%   formulation_comparison\<MODEL>\ and can be viewed with launch_comparison.

    if nargin < 1
        only = '';
    end

    engine_dir = fileparts(mfilename('fullpath'));
    root_dir = fileparts(engine_dir);
    stub_dir = fullfile(engine_dir, 'stub');
    refs_dir = fileparts(root_dir);
    book_dir = fullfile(refs_dir, ...
        'Nikravesh - Planar Multibody Dynamics 2nd Edition');
    bcf_dir = fullfile(book_dir, 'DAP_BC_12_2018', 'DAP_BC 12_2018');
    jcf_dir = fullfile(book_dir, 'DAP-JC_11-2018', 'DAP_JC 11_2018');

    path_orig = path;
    restore_path = onCleanup(@() path(path_orig));
    addpath(engine_dir);
    warn_state = warning('off', 'MATLAB:dispatcher:nameConflict');
    restore_warn = onCleanup(@() warning(warn_state));

    % Columns: output name, BCF model folder, JCF model folder, tf [s], dt [s]
    models = {
        'AA',     'AA',     'AA',  1.0, 0.01
        'Cart_C', 'Cart_C', 'Cart_C', 1.0, 0.01
        'DT',     'DT',     'DT',  5.0, 0.01
        'HN',     'HN',     'HN',  1.0, 0.01
        'MB_A',   'MB_A',   'MBA', 1.0, 0.01
        'MB_B',   'MB_B',   'MBB', 1.0, 0.01
        'MC',     'MC',     'MC',  1.0, 0.01
        'MP_A',   'MP_A',   'MPA', 1.0, 0.01
        'MP_B',   'MP_B',   'MPB', 1.0, 0.01
        'SAT',    'SAT',    'SAT', 3.0, 0.01
        'SP',     'SP',     'SP',  1.0, 0.01
        };

    unpaired = {
        'BCF only: Cart_A, Cart_B, Cart_D, CB, MP_C'
        'JCF only: Rod2, SC, SCd, TP'
        };
    fprintf('Not run (no counterpart in the other formulation):\n');
    fprintf('  %s\n', unpaired{:});

    if ~isempty(only) && ~any(strcmpi(models(:, 1), only))
        error('run_comparisons:unknownModel', ...
            'Model "%s" is not in the paired model table.', only);
    end

    failures = cell(0, 3);
    for k = 1:size(models, 1)
        name = models{k, 1};
        if ~isempty(only) && ~strcmpi(name, only)
            continue
        end

        out_dir = fullfile(root_dir, name);
        if ~isfolder(out_dir)
            mkdir(out_dir);
        end

        forms = {'BCF', bcf_dir, models{k, 2}; ...
                 'JCF', jcf_dir, models{k, 3}};
        for f = 1:size(forms, 1)
            cfg = struct('model', name, 'formulation', forms{f, 1}, ...
                'solver_dir', forms{f, 2}, 'folder', forms{f, 3}, ...
                'tf', models{k, 4}, 'dt', models{k, 5}, 'out_dir', out_dir);

            path_before = path;
            old_dir = cd(cfg.solver_dir);
            addpath(stub_dir, '-begin');
            try
                run_formulation(cfg);
            catch err
                where = '';
                if ~isempty(err.stack)
                    where = sprintf(' [%s, line %d]', ...
                        err.stack(1).name, err.stack(1).line);
                end
                msg = [err.message where];
                failures(end+1, :) = {name, cfg.formulation, msg};
                warning('run_comparisons:failed', ...
                    '%s %s failed: %s', name, cfg.formulation, msg);
            end
            path(path_before);
            cd(old_dir);
        end
    end

    if isempty(failures)
        fprintf('All requested runs completed.\n');
    else
        fprintf('%d run(s) failed:\n', size(failures, 1));
        for r = 1:size(failures, 1)
            fprintf('  %s %s: %s\n', failures{r, :});
        end
    end
end
