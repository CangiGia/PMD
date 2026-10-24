function value = input(prompt, option)
%INPUT Non-interactive replacement for the built-in input function.
%   Answers the prompts of the DAP drivers (dap.m) from the appdata entry
%   'cmp_answers', which is set by run_formulation before dap is called.
%   Only on the MATLAB path while run_comparisons is running.

    answers = getappdata(0, 'cmp_answers');

    if contains(prompt, 'folder', 'IgnoreCase', true)
        value = answers.folder;
    elseif contains(prompt, 'correct', 'IgnoreCase', true)
        value = answers.ic;
    elseif contains(prompt, 'Final time', 'IgnoreCase', true)
        value = answers.tf;
    elseif contains(prompt, 'Reporting time-step', 'IgnoreCase', true)
        value = answers.dt;
    else
        error('input:unhandledPrompt', 'Unhandled prompt: %s', prompt);
    end

    if nargin < 2 || ~strcmp(option, 's')
        value = double(value);
    end
end
