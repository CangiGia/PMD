function add_labels_to_results(only)
%ADD_LABELS_TO_RESULTS Add q_labels and q_kinds to existing result files.
%   add_labels_to_results() appends the physical labels to every
%   <formulation>_results.mat already in formulation_comparison\<MODEL>\.
%   The trajectories are not recomputed, so this is fast. Labels come from
%   coord_labels, the same function that run_formulation uses.
%   add_labels_to_results('MODEL') processes a single model folder.

    if nargin < 1
        only = '';
    end

    engine_dir = fileparts(mfilename('fullpath'));
    root_dir = fileparts(engine_dir);
    refs_dir = fileparts(root_dir);
    book_dir = fullfile(refs_dir, ...
        'Nikravesh - Planar Multibody Dynamics 2nd Edition');
    bcf_dir = fullfile(book_dir, 'DAP_BC_12_2018', 'DAP_BC 12_2018');
    jcf_dir = fullfile(book_dir, 'DAP-JC_11-2018', 'DAP_JC 11_2018');

    % Columns: output folder, BCF model folder, JCF model folder
    models = {
        'AA',     'AA',     'AA'
        'Cart_C', 'Cart_C', 'Cart_C'
        'DT',     'DT',     'DT'
        'HN',     'HN',     'HN'
        'MB_A',   'MB_A',   'MBA'
        'MB_B',   'MB_B',   'MBB'
        'MC',     'MC',     'MC'
        'MP_A',   'MP_A',   'MPA'
        'MP_B',   'MP_B',   'MPB'
        'Rod',    'Rod',    'Rod'
        'SAT',    'SAT',    'SAT'
        'SP',     'SP',     'SP'
        };

    for k = 1:size(models, 1)
        name = models{k, 1};
        if ~isempty(only) && ~strcmpi(name, only)
            continue
        end
        out_dir = fullfile(root_dir, name);
        forms = {'BCF', bcf_dir, models{k, 2}; ...
                 'JCF', jcf_dir, models{k, 3}};
        for f = 1:size(forms, 1)
            res_file = fullfile(out_dir, [forms{f, 1} '_results.mat']);
            if ~isfile(res_file)
                fprintf('%s %s: no results file, skipped.\n', name, forms{f, 1});
                continue
            end
            data = load(res_file, 'q');
            cfg = struct('model', name, 'formulation', forms{f, 1}, ...
                'solver_dir', forms{f, 2}, 'folder', forms{f, 3});
            [q_labels, q_kinds] = coord_labels(size(data.q, 2)); %#ok<ASGLU>
            save(res_file, 'q_labels', 'q_kinds', '-append');
            fprintf('%s %s: %d labels written.\n', name, forms{f, 1}, ...
                numel(q_labels));
        end
    end
end
