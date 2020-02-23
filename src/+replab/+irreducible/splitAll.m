function irreps = splitAll(rep, context)
% Decomposes fully the given representation into subrepresentations
%
% Returns a list of irreducible representations, where trivial subrepresentations
% have been identified
%
% If the given representation is irreducible, will set the `+replab.Rep.isIrreducible`
% attribute and return it.
%
% Args:
%   rep (`+replab.Rep`): Representation to decompose
%   context (`+replab.Context`): Sampling context to use
%
% Returns:
%   cell(1,\*) of `+replab.SubRep`: irreducible subrepresentations
    d = rep.dimension;
    start = replab.rep.fullSubRep(rep);
    todo = {start};
    irreps = cell(1, 0);
    while ~isempty(todo)
        h = todo{1};
        if isequal(h.isIrreducible, true)
            % head of list is irreducible, remove it
            irreps{1,end+1} = h;
            todo = todo(2:end);
        else
            res = replab.irreducible.split(h, context);
            res = cellfun(@(sub) replab.rep.collapseSubRepSubRep(sub), res, 'uniform', 0);
            todo = horzcat(todo(2:end), res);
        end
    end
end
