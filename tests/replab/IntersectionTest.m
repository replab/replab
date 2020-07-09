function test_suite = IntersectionTest()
    disp(['Setting up tests in ', mfilename()]);
    try
        test_functions = localfunctions();
    catch
    end
    initTestSuite;
end

function test_intersection_commutative
    n = 8;
    Sn = replab.S(n);
    f = @() replab.nfg.randomSubgroup(replab.nfg.randomSubgroup(Sn));
    G1 = f();
    G2 = f();
    G3 = f();
    I1 = G1.intersection(G2.intersection(G3));
    I2 = G3.intersection(G1.intersection(G2));
    assert(I1 == I2);
end

function test_intersection_trivial
    S3 = replab.S(3);
    T = S3.trivialSubgroup;
    T.intersection(S3);
end