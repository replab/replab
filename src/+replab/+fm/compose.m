function m = compose(second, first)
% Returns the composition of two finite morphisms
%
% If both morphisms are isomorphisms, returns an isomorphism
%
%
% Args:
%   second (`.FiniteMorphism`): Morphism to apply second
%   first (`.FiniteMorphism`): Morphism to apply first
%
% Returns:
%   `.Morphism`: The composition of the two morphisms
    assert(second.source.isSubgroupOf(first.target));
    if isa(first, 'replab.FiniteIsomorphism') && isa(second, 'replab.FiniteIsomorphism')
        m = replab.fm.IsoComposition(second, first);
    else
        m = replab.fm.Composition(second, first);
    end
end
