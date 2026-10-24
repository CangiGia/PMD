function [labels, kinds] = coord_labels(nq)
%COORD_LABELS Physical name and unit kind of each generalized coordinate.
%   [labels, kinds] = COORD_LABELS(nq) returns one label and one kind per
%   column of q (nq columns, three per body).
%
%   Both formulations store q as (x, y, phi) per body, so the labels are
%   "B<i> x", "B<i> y" and "B<i> phi" in body order. kinds is 'm' for
%   lengths, 'rad' for angles.

    [labels, kinds] = bcf_labels(nq);
end

function [labels, kinds] = bcf_labels(nq)
%BCF_LABELS Labels for the (x, y, phi) body coordinates of BCF.
    comp = {'x', 'y', 'phi'};
    ckind = {'m', 'm', 'rad'};
    labels = cell(1, nq);
    kinds = cell(1, nq);
    for k = 1:nq
        body = ceil(k / 3);
        c = mod(k - 1, 3) + 1;
        labels{k} = sprintf('B%d %s', body, comp{c});
        kinds{k} = ckind{c};
    end
end

