function tuples = cartesian(sets)
% Returns the cartesian product of sets
%
% Args:
%   sets (cell(1,nSets) of cell(1,\*) of elements): Sets to take the cartesian product of
%
%
% Returns:
%   cell(1,\*) of cell(1,nSets): Tuples of elements
    nSets = length(sets);
    nEls = cellfun(@length, sets);
    M = replab.util.MixedRadix(nEls, true, false);
    tuples = cell(1, prod(nEls));
    for i = 1:prod(nEls)
        sub = M.ind2sub(i);
        el = cell(1, nSets);
        for j = 1:nSets
            el{j} = sets{j}{sub(j)};
        end
        tuples{i} = el;
    end
end
