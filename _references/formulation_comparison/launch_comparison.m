function launch_comparison(model)
%LAUNCH_COMPARISON Print statistics and plot the BCF/JCF results of a model.
%   launch_comparison('MODEL') reads MODEL\BCF_results.mat and
%   MODEL\JCF_results.mat (variables q, reac, lambda, t), prints
%   MODEL\stats_BCF.txt and MODEL\stats_JCF.txt to the command window, and
%   plots q, reac and lambda against t for both formulations.
%
%   Model folders are detected automatically: every subfolder of this
%   folder that contains BCF_results.mat or JCF_results.mat is a model.
%   The name is matched case-insensitively.
%
%   launch_comparison() with no argument, or with an unknown name, prints
%   the list of available models and stops with an error.
%
%   Figures:
%     Figure 1 - q versus t. One row per formulation (BCF on top, JCF
%                below). Both formulations report (x, y, phi) per body,
%                so the panels line up. Columns: positions [m] and angles
%                [deg]; a third column appears only when some coordinate
%                has an unknown unit. A legend is omitted above 24 lines.
%     Figure 2 - reac versus t: constraint force on each body. One row per
%                formulation. Columns: forces [N] (Fx, Fy) and moments
%                [N m] (Mz). BCF is D'*Lambda, JCF is M*c_dd - h_a.
%     Figure 3 - lambda versus t (BCF on top, JCF below). Raw multipliers.
%                A panel reads "no constraints" when the multiplier matrix
%                is empty, which is the case for an open JCF chain.


    root = fileparts(mfilename('fullpath'));
    models = find_models(root);

    if nargin < 1 || isempty(model)
        print_models(models);
        error('launch_comparison:noModel', ...
            'Specify a model name, for example launch_comparison(''%s'').', ...
            first_or_placeholder(models));
    end

    idx = find(strcmpi(models, model), 1);
    if isempty(idx)
        print_models(models);
        error('launch_comparison:unknownModel', ...
            'Model "%s" was not found in %s.', model, root);
    end
    name = models{idx};
    model_dir = fullfile(root, name);

    forms = {'BCF', 'JCF'};
    data = struct();
    for f = 1:numel(forms)
        fm = forms{f};
        fprintf('==================== %s - %s ====================\n', ...
            name, fm);
        stats_file = fullfile(model_dir, ['stats_' fm '.txt']);
        if isfile(stats_file)
            fprintf('%s\n', fileread(stats_file));
        else
            fprintf('Statistics file not found: %s\n\n', stats_file);
        end

        res_file = fullfile(model_dir, [fm '_results.mat']);
        if isfile(res_file)
            data.(fm) = load(res_file);
        end
    end

    fig_q = figure('Name', [name ' - generalized coordinates q'], ...
        'NumberTitle', 'off', 'Position', [50 50 1400 800]);
    col_titles = {'positions', 'angles', 'unit unknown'};
    if has_unknown_q(data, forms)
        ncols = numel(col_titles);
    else
        ncols = 2;
    end
    for f = 1:numel(forms)
        fm = forms{f};
        for c = 1:ncols
            subplot(2, ncols, (f - 1) * ncols + c, 'Parent', fig_q);
            ttl = sprintf('%s: %s', fm, col_titles{c});
            if isfield(data, fm)
                draw_q(data.(fm), ttl, c);
            else
                draw_message(sprintf('%s results not found', fm), ttl);
            end
        end
    end

    fig_r = figure('Name', [name ' - constraint forces on bodies'], ...
        'NumberTitle', 'off', 'Position', [100 100 1200 800]);
    reac_titles = {'forces [N]', 'moments [N m]'};
    for f = 1:numel(forms)
        fm = forms{f};
        for c = 1:2
            subplot(2, 2, (f - 1) * 2 + c, 'Parent', fig_r);
            ttl = sprintf('%s: %s', fm, reac_titles{c});
            if isfield(data, fm)
                draw_reac(data.(fm), ttl, c);
            else
                draw_message(sprintf('%s results not found', fm), ttl);
            end
        end
    end

    fig_l = figure('Name', [name ' - Lagrange multipliers lambda'], ...
        'NumberTitle', 'off');
    for f = 1:numel(forms)
        fm = forms{f};
        subplot(2, 1, f, 'Parent', fig_l);
        if isfield(data, fm)
            draw_panel(data.(fm).t, data.(fm).lambda, ...
                sprintf('%s: lambda(t)', fm), 'lambda', 'lambda');
        else
            draw_message(sprintf('%s results not found', fm), ...
                sprintf('%s: lambda(t)', fm));
        end
    end
end

function tf = has_unknown_q(data, forms)
%HAS_UNKNOWN_Q True if any loaded formulation has a q column of unknown unit.
%   Files without q_kinds count as unknown, so the third column is kept.
    tf = false;
    for f = 1:numel(forms)
        fm = forms{f};
        if ~isfield(data, fm)
            continue
        end
        res = data.(fm);
        if ~isfield(res, 'q_kinds') || any(strcmp(res.q_kinds, '?'))
            tf = true;
            return
        end
    end
end

function names = find_models(root)
%FIND_MODELS Return the names of the model subfolders of ROOT.
    entries = dir(root);
    names = {};
    for k = 1:numel(entries)
        e = entries(k);
        if ~e.isdir || startsWith(e.name, '.') || strcmp(e.name, 'engine')
            continue
        end
        model_dir = fullfile(root, e.name);
        if isfile(fullfile(model_dir, 'BCF_results.mat')) || ...
                isfile(fullfile(model_dir, 'JCF_results.mat'))
            names{end+1} = e.name; %#ok<AGROW>
        end
    end
end

function print_models(models)
%PRINT_MODELS Print the detected model names.
    if isempty(models)
        fprintf('No model folders with results were found.\n');
    else
        fprintf('Available models:\n');
        fprintf('  %s\n', models{:});
    end
end

function name = first_or_placeholder(models)
%FIRST_OR_PLACEHOLDER Return the first model name, or a generic placeholder.
    if isempty(models)
        name = 'MODEL';
    else
        name = models{1};
    end
end

function draw_panel(t, y, ttl, ylab, prefix)
%DRAW_PANEL Plot the columns of Y against T, or show "no constraints".
    if isempty(y)
        draw_message('no constraints', ttl);
        return
    end
    n = size(y, 2);
    set(gca, 'ColorOrder', turbo(n));
    plot(t(:), y);
    if n <= 24
        labels = arrayfun(@(i) sprintf('%s_%d', prefix, i), ...
            1:n, 'UniformOutput', false);
        legend(labels, 'Location', 'eastoutside', 'Interpreter', 'none');
    end
    title(ttl);
    xlabel('t [s]');
    ylabel(ylab);
    grid on;
end

function draw_q(res, ttl, group)
%DRAW_Q Plot one unit group of q against t with physical legend names.
%   GROUP is 1 for lengths [m], 2 for angles (plotted in deg) and 3 for
%   coordinates with an unknown unit. Files written before the labels
%   were added have no q_labels or q_kinds, so every column falls in
%   group 3 with a generic name.
    t = res.t(:);
    q = res.q;
    n = size(q, 2);
    if isfield(res, 'q_labels') && isfield(res, 'q_kinds')
        labels = res.q_labels;
        kinds = res.q_kinds;
    else
        labels = arrayfun(@(k) sprintf('q_%d', k), 1:n, ...
            'UniformOutput', false);
        kinds = repmat({'?'}, 1, n);
    end

    switch group
        case 1
            draw_columns(t, q, labels, kinds, ttl, 'm', ...
                'position [m]', 1);
        case 2
            draw_columns(t, q, labels, kinds, ttl, 'rad', ...
                'angle [deg]', 180 / pi);
        otherwise
            draw_columns(t, q, labels, kinds, ttl, '?', ...
                'q [raw]', 1);
    end
end

function draw_reac(res, ttl, group)
%DRAW_REAC Plot the constraint force (group 1) or moment (group 2).
%   Both formulations save reac with reac_labels and reac_kinds. Files
%   written before that change have neither, so the panel says so.
    if ~isfield(res, 'reac') || ~isfield(res, 'reac_kinds')
        draw_message('reac not saved, rerun run_comparisons', ttl);
        return
    end
    t = res.t(:);
    if group == 1
        draw_columns(t, res.reac, res.reac_labels, res.reac_kinds, ...
            ttl, 'N', 'force [N]', 1);
    else
        draw_columns(t, res.reac, res.reac_labels, res.reac_kinds, ...
            ttl, 'Nm', 'moment [N m]', 1);
    end
end

function draw_columns(t, Y, labels, kinds, ttl, want, ylab, scale)
%DRAW_COLUMNS Plot the columns of Y whose kind is WANT against t.
%   The columns are multiplied by SCALE. The legend is omitted above 24
%   lines, and the title then gives the number of lines instead.
    idx = find(strcmp(kinds, want));
    if isempty(idx)
        draw_message('none', ttl);
        return
    end

    y = Y(:, idx) * scale;
    set(gca, 'ColorOrder', turbo(numel(idx)));
    plot(t, y);
    if numel(idx) <= 24
        legend(labels(idx), 'Location', 'eastoutside', ...
            'Interpreter', 'none');
        title(ttl);
    else
        title(sprintf('%s (%d lines, legend omitted)', ttl, numel(idx)));
    end
    xlabel('t [s]');
    ylabel(ylab);
    grid on;
end


function draw_message(msg, ttl)
%DRAW_MESSAGE Show a centred text message in an empty panel.
    axis off;
    text(0.5, 0.5, msg, 'Units', 'normalized', ...
        'HorizontalAlignment', 'center');
    title(ttl);
end
