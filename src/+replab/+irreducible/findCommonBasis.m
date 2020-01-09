function W = findCommonBasis(rep, samples, sub1, sub2)
% Finds the change of basis that relates two equivalent subrepresentations
%
% Args:
%   rep (replab.Rep): Representation decomposed
%   samples (replab.irreducible.OnDemandSamples): Lazy evaluation of various samples for ``rep``
%   sub1 (replab.SubRep): Irreducible subrepresentation of rep
%   sub2 (replab.SubRep): Irreducible subrepresentation of rep equivalent to ``sub1``
%
% Returns ``W`` such that ``W * sub2.image(g) * W' = sub1.image(g)``
    assert(isa(sub1, 'replab.SubRep'));
    assert(isa(sub2, 'replab.SubRep'));
    assert(sub1.parent == rep);
    assert(sub2.parent == rep);
    assert(sub1.isKnownIrreducible);
    assert(sub2.isKnownIrreducible);
    d = sub1.dimension;
    assert(sub2.dimension == d);
    c = samples.commutantSample(2);
    W = sub1.U*c*sub2.U';
    % correct the normalization
    W = W * sqrt(d/real(trace(W*W')));
end
