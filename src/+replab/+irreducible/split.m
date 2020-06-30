function sub = split(rep, context)
% Decomposes the given representation into subrepresentations
%
% The returned list of subrepresentations contains instances of `+replab.SubRep`, whose
% `~+replab.SubRep.parent` property must be equal to the given representation ``rep``.
%
% If the returned list contains a single subrepresentation, it needs to have a trivial
% basis and embedding map equal to the identity matrix, and that subrepresentation needs
% to be more informative than the given representation by either having `~+replab.Rep.isIrreducible`
% or `~+replab.Rep.trivialDimension` filled up.
%
% If the return list contains more than one subrepresentation, no such restriction applies.
%
% See the default implementations of that method in `+replab.dispatchDefaults`
%
% Args:
%   rep (`+replab.Rep`): Representation to decompose, must be not known to be irreducible
%   context (`+replab.Context`): A context in which to cache samples
%
% Returns:
%   cell(1,\*) of `+replab.SubRep`: A cell array of subrepresentations
    assert(isa(rep, 'replab.Rep'));
    assert(isa(context, 'replab.Context'));
    if ~replab.dispatch('exists', 'replab.irreducible.split')
        replab.dispatch('register', 'replab.irreducible.split', 'splitTrivialGroup', 500, ...
                        @(r, c) replab.irreducible.splitTrivialGroup(r, c));
        replab.dispatch('register', 'replab.irreducible.split', 'splitBlocks', 400, ...
                        @(r, c) replab.irreducible.splitBlocks(r, c));
        replab.dispatch('register', 'replab.irreducible.split', 'splitAllOnes', 300, ...
                        @(r, c) replab.irreducible.splitAllOnes(r, c));
        replab.dispatch('register', 'replab.irreducible.split', 'splitTrivialComponent', 200, ...
                        @(r, c) replab.irreducible.splitTrivialComponent(r, c));
        replab.dispatch('register', 'replab.irreducible.split', 'splitUsingCommutant', 100, ...
                        @(r, c) replab.irreducible.splitUsingCommutant(r, c));
    end
    sub = replab.dispatch('call', 'replab.irreducible.split', rep, context);
    if isa(sub, 'replab.DispatchNext')
        return
    end
    for i = 1:length(sub)
        assert(isa(sub{i}, 'replab.SubRep'));
        assert(sub{i}.parent == rep);
    end
end
