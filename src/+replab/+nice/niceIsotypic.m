function iso1 = niceIsotypic(iso)
% Attempts to make the bases of an isotypic component nicer
%
% Args:
%   iso (`+replab.Isotypic`): Isotypic component to work on
%
% Returns:
%   `+replab.Isotypic` or `+replab.DispatchNext`: Isotypic component with nicer basis if found
    assert(isa(iso, 'replab.Isotypic'));
    iso1 = replab.dispatch('call', 'replab.nice.niceIsotypic', iso);
    if isa(iso1, 'replab.DispatchNext')
        return
    end
    assert(isa(iso1, 'replab.SubRep'));
    assert(iso1.parent == iso.parent);
end
