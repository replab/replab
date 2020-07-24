function test_suite = ConjugacyClassesTest()
    disp(['Setting up tests in ', mfilename()]);
    try
        test_functions = localfunctions();
    catch
    end
    initTestSuite;
end


function test_conjugacyClassRepresentative
    if ReplabTestParameters.onlyFastTests
        n = 1;
        d = 4;
    else
        n = 10;
        d = 15;
    end
    S = replab.S(d);
    for i = 1:n
        G = S.randomProperSubgroup(3);
        g1 = G.sample;
        s = G.sample;
        g2 = G.leftConjugate(s, g1);
        r1 = G.conjugacyClass(g1).representative;
        r2 = G.conjugacyClass(g2).representative;
        assertEqual(r1, r2);
    end
end
