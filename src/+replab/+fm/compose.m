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
    if isa(first, 'replab.fm.Identity')
        m = second;
        return
    end
    if isa(second, 'replab.fm.Identity')
        m = first;
        return
    end
    areIsomorphisms = isa(first, 'replab.Isomorphism') && isa(second, 'replab.Isomorphism');
    areFinite = isa(first, 'replab.FiniteMorphism') && isa(second, 'replab.FiniteMorphism');
    if areIsomorphisms
        if areFinite
            m = replab.fm.FiniteIsomorphismComposition(second, first);
        else
            m = replab.fm.IsomorphismComposition(second, first);
        end
    else
        if areFinite
            m = replab.fm.FiniteComposition(second, first);
        else
            m = replab.fm.Composition(second, first);
        end
    end
end
