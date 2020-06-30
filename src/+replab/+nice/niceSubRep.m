function sub1 = niceSubRep(sub)
% Attempts to make the basis of a subrepresentation nicer
%
% Args:
%   rep (`+replab.SubRep`): Subrepresentation to work on
%
% Returns:
%   `+replab.SubRep` or `+replab.DispatchNext`: Subrepresentation with nicer basis if found
    assert(isa(sub, 'replab.SubRep'));
    if ~replab.dispatch('exists', 'replab.nice.niceSubRep')
        replab.dispatch('register', 'replab.nice.niceSubRep', 'recoverReal', 10, ...
                        @(sub) replab.nice.niceSubRepRecoverReal(sub));
        replab.dispatch('register', 'replab.nice.niceSubRep', 'recoverInteger', 5, ...
                        @(sub) replab.nice.niceSubRepRecoverInteger(sub));
        replab.dispatch('register', 'replab.nice.niceSubRep', 'original', 0, ...
                        @(sub) replab.nice.niceSubRepOriginal(sub));
    end
    sub1 = replab.dispatch('call', 'replab.nice.niceSubRep', sub);
    if isa(sub1, 'replab.DispatchNext')
        return
    end
    assert(isa(sub1, 'replab.SubRep'));
    assert(sub1.parent == sub.parent);
end
