function W = findCommonBasis(rep, sub1, sub2, c)
% Finds the change of basis that relates two copies of an irrep
%
%      rep: Representation decomposed
%     sub1: First copy of an irreducible subrepresentation of rep
%           whose division algebra structure is canonical
%           (instance of replab.Irrep)
%     sub2: Second copy of an irreducible subrepresentation of rep
%           (sub1 and sub2 need to be equivalent)
%        c: (optional) generic sample from rep.commutant
%
% Returns W such that W * sub2.image(g) * W' = sub1.image(g)
    assert(isa(sub1, 'replab.Irrep'));
    assert(isa(sub2, 'replab.SubRep'));
    if nargin < 4
        c = rep.commutant.sample;
    end
    d = sub1.dimension;
    assert(sub2.dimension == d);
    W = sub2.U*c*sub1.U';
    W = W * sqrt(sub1.dimension/real(trace(W*W')));
    % correction needed, TODO: why?
    W = W';
end
