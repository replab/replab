function S = SignedSymmetricGroup(n)
    T = replab.signed.FiniteGroupType.make(n);
    S = T.isomorphism.source;
end
