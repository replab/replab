function test_suite = SubRepTest()
    disp(['Setting up tests in ', mfilename()]);
    try
        test_functions = localfunctions();
    catch
    end
    initTestSuite;
    F = [ 2  0
          -1  1
          -1 -1];
    G = [2 -1 -1
         0  3 -3]/6;
    S3 = replab.S(3);
    rep = S3.naturalRep.subRep(F, 'projection', G);
    test_suite = rep.laws.addTestCases(test_suite);
end

function test_quaterion_refine
    Q = replab.QuaternionGroup();
    S3 = replab.S(3);
    W = S3.wreathProduct(Q);
    rep = W.primitiveRep(Q.naturalRep);
    irreps = rep.split;
    sub = irreps{1};
    subN = sub.withNoise(0.001, 0.001);
    gen = replab.rep.GenSubRep.fromQuaternionTypeSubRep(subN);

    tol = replab.rep.Tolerances;

    gen1 = replab.rep.refine_nonUnitary_largeScale(gen, 5, tol, [], []);
    sub1 = gen1.toSubRep;

    gen2 = replab.rep.refine_unitary_largeScale(gen, 5, tol, []);
    sub2 = gen2.toSubRep;

    gen3 = replab.rep.refine_unitary_mediumScale(gen, 5, tol, []);
    sub3 = gen3.toSubRep;

    gen4 = replab.rep.refine_nonUnitary_mediumScale(gen, 3, tol, [], []);
    sub4 = gen4.toSubRep;

    gen1.check
    gen2.check
    gen3.check
    gen4.check

    sub1.check
    sub2.check
    sub3.check
    sub4.check

    assert(strcmp(sub1.divisionAlgebraName, 'quaternion.rep'));
    assert(strcmp(sub2.divisionAlgebraName, 'quaternion.rep'));

end

function test_complex_refine
    g1 = [6 3 4 2 5 7 1];
    g2 = [6 4 1 7 2 5 3];

    group = replab.PermutationGroup.of(g1, g2);
    rep = group.naturalRep;
    dec = rep.decomposition;
    sub = dec.irrep(2);
    subN = sub.withNoise(0.01);
    gen = replab.rep.GenSubRep.fromComplexTypeSubRep(subN);

    tol = replab.rep.Tolerances;

    gen1 = replab.rep.refine_nonUnitary_largeScale(gen, 5, tol, [], []);
    sub1 = gen1.toSubRep;

    gen2 = replab.rep.refine_unitary_largeScale(gen, 5, tol, []);
    sub2 = gen2.toSubRep;

    gen3 = replab.rep.refine_unitary_mediumScale(gen, 5, tol, []);
    sub3 = gen3.toSubRep;

    gen4 = replab.rep.refine_nonUnitary_mediumScale(gen, 5, tol, [], []);
    sub4 = gen4.toSubRep;

    assert(sub1.projectorErrorBound < 1e-14);
    assert(sub2.projectorErrorBound < 1e-14);
    assert(sub3.projectorErrorBound < 1e-14);
    assert(sub4.projectorErrorBound < 1e-14);
    gen1.check
    gen2.check
    gen3.check
    gen4.check

    sub1.check
    sub2.check
    sub3.check
    sub4.check

    assert(strcmp(sub1.divisionAlgebraName, 'complex'));
    assert(strcmp(sub2.divisionAlgebraName, 'complex'));
    assert(strcmp(sub3.divisionAlgebraName, 'complex'));
    assert(strcmp(sub4.divisionAlgebraName, 'complex'));
end

function test_real_refine
    n = 5;
    G = replab.S(n);
    nat = G.naturalRep;
    triv = nat.subRep(ones(n, 1));
    std = nat.maschke(triv);
    pow4 = nat.tensorPower(4);
    sub = pow4.subRep(kron(kron(triv.injection, std.injection), kron(triv.injection, std.injection)));
    subNoisy = sub.withNoise(0.01);
    % the error bound is horrific
    % fprintf('Error bound with noise: %e\n', subNoisy.errorBound);
    sub1 = subNoisy.refine;
    sub2 = subNoisy.refine('largeScale', true, 'nSamples', 5);
    assert(sub1.projectorErrorBound < 1e-14);
    assert(sub2.projectorErrorBound < 1e-14);
    sub1.check
    sub2.check
end
