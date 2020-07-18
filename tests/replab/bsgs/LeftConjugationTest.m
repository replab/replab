function test_suite = LeftConjugationTest()
    disp(['Setting up tests in ', mfilename()]);
    try
        test_functions = localfunctions();
    catch
    end
    initTestSuite;
end

function test_findLeftConjugations
    if ReplabTestParameters.onlyFastTests
        n = 1;
        d = 3;
    else
        n = 10;
        d = 10;
    end
    for i = 1:n
        S = replab.S(d);
        G = S.randomProperSubgroup(2);
        s = G.sample;
        g = G.sample;
        t = G.leftConjugate(g, s);
        L = G.findLeftConjugations(s, t);
        l = L.sample;
        G.assertEqv(G.leftConjugate(l, s), t);
    end
end
