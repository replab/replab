Q = replab.QuaternionGroup();
S3 = replab.S(3);
W = S3.wreathProduct(Q);
rep = W.primitiveRep(Q.naturalRep);
irreps = rep.split;
sub = irreps{1};
subN = sub.withNoise(0.1);
gen = replab.GenSubRep.fromQuaternionTypeSubRep(subN);

gen1 = replab.rep.refine_nonUnitary_genSubRep_largeScale(gen, 20, 5, 100, [], []);
sub1 = gen1.toSubRep;

gen2 = replab.rep.refine_unitary_genSubRep_largeScale(gen, 20, 5, 100, []);
sub2 = gen2.toSubRep;

gen3 = replab.rep.refine_unitary_genSubRep_mediumScale(gen, 3, 10, []);
sub3 = gen3.toSubRep;

gen4 = replab.rep.refine_nonUnitary_genSubRep_mediumScale(gen, 3, 10, []);
sub4 = gen3.toSubRep;
