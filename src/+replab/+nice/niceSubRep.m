function sub1 = niceSubRep(sub)
% Attempts to make the basis of a subrepresentation nicer
%
% Args:
%   rep (`+replab.SubRep`): Subrepresentation to work on
%
% Returns:
%   `+replab.SubRep` or `+replab.DispatchNext`: Subrepresentation with nicer basis if found
    assert(isa(sub, 'replab.SubRep'));
    sub1 = replab.dispatch('call', 'replab.nice.niceSubRep', sub);
    if isa(sub1, 'replab.DispatchNext')
        return
    end
    assert(isa(sub1, 'replab.SubRep'));
    assert(sub1.parent == sub.parent);
end
