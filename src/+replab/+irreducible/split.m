function sub = split(rep, context)
% Decomposes the given representation into subrepresentations
%
% Must return either:
% - A cell array of two or more subrepresentations identified in ``rep``, they must
%   be of type `replab.SubRep` and have ``rep`` as parent,
% - A cell array of a single representation, where ``sub{1} == rep``, but at least
%   one of ``rep.isIrreducible`` or ``rep.trivialDimension`` have been computed
%   whereas they were empty before.
%
% Note: see the default methods applied in `+replab.dispatchDefaults`
%
% Args:
%   rep (`+replab.Rep`): Representation to decompose, must be not known to be irreducible
%   context (`+replab.Context`): A context in which to cache samples
%
% Returns:
%   cell(1,*) of `+replab.Rep` or `+replab.SubRep`: A cell array of subrepresentations
    assert(isa(rep, 'replab.Rep'));
    assert(isa(context, 'replab.Context'));
    sub = replab.dispatch('call', 'replab.irreducible.split', rep, context);
end
