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
    if isa(first, 'replab.mrp.Identity')
        m = second;
        return
    end
    if isa(second, 'replab.mrp.Identity')
        m = first;
        return
    end
    if isa(second, 'replab.Rep')
        m = second.compose(first);
        return
    end
    areIsomorphisms = isa(first, 'replab.Isomorphism') && isa(second, 'replab.Isomorphism');
    areFinite = isa(first, 'replab.FiniteMorphism') && isa(second, 'replab.FiniteMorphism');
    if areIsomorphisms
        if areFinite
            m = replab.mrp.FiniteIsomorphismComposition(second, first);
        else
            m = replab.mrp.IsomorphismComposition(second, first);
        end
    else
        if areFinite
            m = replab.mrp.FiniteComposition(second, first);
        else
            m = replab.mrp.Composition(second, first);
        end
    end
end
