function m = compose(second, first)
% Returns the composition of two finite morphisms
%
% Returns the most general type afforded by the arguments.
%
% Args:
%   second (`+replab.Morphism`): Morphism to apply second
%   first (`+replab.Morphism`): Morphism to apply first
%
% Returns:
%   `+replab.Morphism`: The composition of the two morphisms
    assert(second.source.isSubgroupOf(first.target));
    if isa(first, 'replab.FiniteIsomorphism') && isa(second, 'replab.FiniteIsomorphism')
        m = replab.fm.IsoComposition(second, first);
    elseif isa(first, 'replab.FiniteMorphism') && isa(second, 'replab.FiniteMorphism')
        m = replab.fm.FiniteComposition(second, first);
    else
        m = replab.fm.Composition(second, first);
    end
end
