function sub1 = refineSubRep(sub)
% Implementation of `+replab.SubRep.refine`
    d = sub.dimension;
    P = full(sub.B_internal*sub.E_internal);
    P = sub.parent.commutant.project(P);
    if isequal(sub.isUnitary, true) && isequal(sub.parent.isUnitary, true)
        [U D] = replab.numerical.sortedEig((P + P')/2, 'descend', false);
        basis = U(:,1:d);
        sub1 = sub.parent.subRep(basis, basis');
    else
        [U,S,V] = svd(P);
        V = V';
        %P = U*S*V
        embedding = S(1:d, 1:d) * V(1:d,:);
        basis = U(:,1:d);
        sub1 = sub.parent.subRep(basis, embedding);
    end
end
