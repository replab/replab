function test_suite = HarmonizeTest()
    disp(['Setting up tests in ', mfilename()]);
    try
        test_functions = localfunctions();
    catch
    end
    initTestSuite;
end

function test_harmonize_quaternion
    Q = replab.QuaternionGroup();
    S3 = replab.S(3);
    W = S3.wreathProduct(Q);
    rep = W.primitiveRep(Q.naturalRep);
    irreps = rep.split;
    dims = cellfun(@(i) i.dimension, irreps);
    inds = find(dims == 16);
    s1 = irreps{inds(1)};
    s2 = irreps{inds(2)};
    g1 = replab.rep.GenSubRep.fromSubRep(s1);
    g2 = replab.rep.GenSubRep.fromSubRep(s2);
    tol = replab.rep.Tolerances;
    g2a = replab.rep.harmonize_nonUnitary_largeScale(g2, g1, 5, tol, s1.injection, s1.projection);
    g2b = replab.rep.harmonize_nonUnitary_largeScale(g2, g1, 5, tol, [], []);
    g2c = replab.rep.harmonize_unitary_largeScale(g2, g1, 5, tol, []);
    w = W.sample;
    img1 = g1.image(w);
    img2a = g2a.image(w);
    img2b = g2b.image(w);
    img2c = g2c.image(w);
    assert(norm(img1 - img2a, 'fro') < 1e-13);
    assert(norm(img1 - img2b, 'fro') < 1e-13);
    assert(norm(img1 - img2c, 'fro') < 1e-13);
end