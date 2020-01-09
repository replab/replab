function sub = split(rep, samples, sub)
% Decomposes the given representation into subrepresentations
%
% Note: see the default methods applied in `+replab.dispatchDefaults`
%
% Args:
%   rep (replab.Rep): Unitary representation to decompose
%   samples (replab.irreducible.OnDemandSamples): Lazy evaluation of various samples
%   sub (replab.SubRep): Subrepresentation of `rep` in which to extract irreducible subrepresentations
%
% Returns:
%   row cell array of replab.SubRep: A cell array of irreducible subrepresentations
%                                    with trivial representations identified with the label '1'
    assert(isa(rep, 'replab.Rep'));
    assert(isa(samples, 'replab.irreducible.OnDemandSamples'));
    assert(isa(sub, 'replab.SubRep'));
    assert(isequal(rep.isUnitary, true));
    assert(sub.parent == rep);
    replab.irreducible.tell('dispatch split dimension %d', sub.dimension);
    sub = replab.dispatch('call', 'replab.irreducible.split', rep, samples, sub);
end
