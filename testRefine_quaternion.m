Q = replab.QuaternionGroup();
S3 = replab.S(3);
W = S3.wreathProduct(Q);
rep = W.primitiveRep(Q.naturalRep);
irreps = rep.split;
sub = irreps{1};
subN = sub.withNoise(0.001, 0.001);
gen = replab.rep.GenSubRep.fromQuaternionTypeSubRep(subN);

gen1 = replab.rep.refine_nonUnitary_largeScale(gen, 20, 5, 100, [], []);
sub1 = gen1.toSubRep;

gen2 = replab.rep.refine_unitary_largeScale(gen, 20, 5, 100, []);
sub2 = gen2.toSubRep;

gen3 = replab.rep.refine_unitary_mediumScale(gen, 3, 10, []);
sub3 = gen3.toSubRep;

gen4 = replab.rep.refine_nonUnitary_mediumScale(gen, 3, 10, []);
sub4 = gen4.toSubRep;

gen1.check
gen2.check
gen3.check
gen4.check

sub1.check
sub2.check
sub3.check
sub4.check
