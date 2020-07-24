for i = 1:100
    i
    G = replab.S(10).randomProperSubgroup(3);
    g1 = G.sample;
    g2 = G.leftConjugate(G.sample, g1);
    r1 = replab.bsgs.ConjugacyClasses.representative(G, g1);
    r2 = replab.bsgs.ConjugacyClasses.representative(G, g2);
    assert(all(r1 == r2));
end
