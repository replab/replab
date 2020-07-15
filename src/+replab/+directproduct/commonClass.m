function className = commonClass(groups)
% Returns the lowest denominator class type for the factor groups
    assert(all(cellfun(@(g) isa(g, 'replab.CompactGroup'), groups)));
    allFinite = all(cellfun(@(g) isa(g, 'replab.FiniteGroup'), groups));
    if allFinite
        className = 'replab.FiniteGroup';
    else
        className = 'replab.CompactGroup';
    end
end
