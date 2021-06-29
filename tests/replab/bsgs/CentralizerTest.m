function test_suite = CentralizerTest()
    disp(['Setting up tests in ', mfilename()]);
    try
        test_functions = localfunctions();
    catch
    end
    initTestSuite;
end


function test_resultUnderConjugation
    if ReplabTestParameters.onlyFastTests
        n = 1;
        d = 3;
    else
        n = 10;
        d = 8;
    end
    S = replab.S(d);
    for i = 1:n
        G = S.randomProperSubgroup(2);
        H = S.randomProperSubgroup(2);
        s = G.sample;
        b0 = G.sample;
        t = S.leftConjugate(b0, s);
        sCentralizer = G.centralizer(s);
        tCentralizer = G.centralizer(t);
        assert(sCentralizer.leftConjugateGroup(b0) == tCentralizer);
    end
end
