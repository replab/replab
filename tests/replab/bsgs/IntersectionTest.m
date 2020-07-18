function test_suite = IntersectionTest()
    disp(['Setting up tests in ', mfilename()]);
    try
        test_functions = localfunctions();
    catch
    end
    initTestSuite;
end

function test_resultUnderConjugation
    n = 10;
    d = 8;
    S = replab.S(d);
    for i = 1:n
        G = S.randomProperSubgroup(2);
        H = S.randomProperSubgroup(2);
        b0 = S.sample;
        GHa = G.intersection(H);
        G1 = G.leftConjugateGroup(b0);
        H1 = H.leftConjugateGroup(b0);
        G1H1 = G1.intersection(H1);
        GHb = G1H1.leftConjugateGroup(S.inverse(b0));
        assert(GHb == GHa);
    end
end
