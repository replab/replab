function W = findCommonBasis(rep, sub1, sub2, context)
% Finds the change of basis that relates two equivalent subrepresentations
%
% We must have ``sub1.parent == rep`` and ``sub2.parent == rep``.
%
% Args:
%   rep (`+replab.Rep`): Representation decomposed
%   sub1 (replab.SubRep): Irreducible subrepresentation of rep
%   sub2 (replab.SubRep): Irreducible subrepresentation of rep
%
% Returns:
%   double(\*,\*): ``W`` such that ``W * sub2.image(g) * W' = sub1.image(g)``
    assert(isa(sub1, 'replab.SubRep'));
    assert(isa(sub2, 'replab.SubRep'));
    assert(sub1.parent == rep);
    assert(sub2.parent == rep);
    assert(isequal(sub1.isIrreducible, true));
    assert(isequal(sub2.isIrreducible, true));
    d = sub1.dimension;
    assert(sub2.dimension == d);
    c = samples.commutantSample(2);
    W = sub1.U*c*sub2.U';
    % correct the normalization
    W = W * sqrt(d/real(trace(W*W')));
end
