function [labels, kinds] = reaction_labels(nr)
%REACTION_LABELS Name and unit kind of each constraint-force column.
%   [labels, kinds] = REACTION_LABELS(nr) returns one label and one kind
%   per column of reac (nr columns, three per body).
%
%   The columns are the constraint force along x and y, and the moment
%   about z, acting on each body: "B<i> Fx" [N], "B<i> Fy" [N] and
%   "B<i> Mz" [N*m]. Both formulations fill reac with the same quantity,
%   so the labels apply to BCF and JCF alike.

    comp = {'Fx', 'Fy', 'Mz'};
    ckind = {'N', 'N', 'Nm'};
    labels = cell(1, nr);
    kinds = cell(1, nr);
    for k = 1:nr
        body = ceil(k / 3);
        c = mod(k - 1, 3) + 1;
        labels{k} = sprintf('B%d %s', body, comp{c});
        kinds{k} = ckind{c};
    end
end
