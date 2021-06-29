function iso1 = niceIsotypic(iso)
% Attempts to make the bases of an isotypic component nicer
%
% Args:
%   iso (`+replab.Isotypic`): Isotypic component to work on
%
% Returns:
%   `+replab.Isotypic` or `+replab.DispatchNext`: Isotypic component with nicer basis if found
    assert(isa(iso, 'replab.Isotypic'));
    if ~replab.dispatch('exists', 'replab.nice.niceIsotypic')
        replab.dispatch('register', 'replab.nice.niceIsotypic', 'trivial', 10, ...
                        @(iso) replab.nice.niceIsotypicTrivial(iso));
        replab.dispatch('register', 'replab.nice.niceIsotypic', 'singleMultiplicity', 5, ...
                        @(iso) replab.nice.niceIsotypicSingleMultiplicity(iso));
        replab.dispatch('register', 'replab.nice.niceIsotypic', 'original', 0, ...
                        @(iso) replab.nice.niceIsotypicOriginal(iso));
    end
    iso1 = replab.dispatch('call', 'replab.nice.niceIsotypic', iso);
    if isa(iso1, 'replab.DispatchNext')
        return
    end
    assert(isa(iso1, 'replab.SubRep'));
    assert(iso1.parent == iso.parent);
end
