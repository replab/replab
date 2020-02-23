function sub = splitTrivialGroup(rep, context)
% If a representation is defined on the trivial group, it is trivial, split accordingly
    if ~isa(rep.group, 'replab.FiniteGroup') || ~rep.group.isTrivial
        sub = replab.DispatchNext('Not defined on a trivial group');
        return
    end
    d = rep.dimension;
    sub = cell(1, d);
    for i = 1:rep.dimension
        B_internal = sparse(i, 1, 1, d, 1);
        E_internal = B_internal';
        s = replab.SubRep(rep, B_internal, E_internal);
        s.trivialDimension = 1;
        s.frobeniusSchurIndicator = 1;
        s.isIrreducible = true;
        sub{1, i} = s;
    end
end
