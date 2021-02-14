S9 = replab.S(9);
group = S9.subgroup({[7 4 1 9 5 2 6 8 3] [7 3 4 2 5 6 9 8 1]});
rep = group.naturalRep;
dec = rep.decomposition;
sub = dec.irrep(2);
subN = sub.withNoise(0.01);
gen = replab.GenSubRep.fromComplexTypeSubRep(subN);

gen1 = replab.rep.refine_nonUnitary_genSubRep_largeScale(gen, 20, 5, 100, [], []);
sub1 = gen1.toSubRep;

gen2 = replab.rep.refine_unitary_genSubRep_largeScale(gen, 20, 5, 100, []);
sub2 = gen2.toSubRep;

gen3 = replab.rep.refine_unitary_genSubRep_mediumScale(gen, 3, 10, []);
sub3 = gen3.toSubRep;

gen4 = replab.rep.refine_nonUnitary_genSubRep_mediumScale(gen, 3, 10, []);
sub4 = gen3.toSubRep;

%sub1 = sub.withUpdatedMaps(I, P, 'divisionAlgebraName', 'complex');
