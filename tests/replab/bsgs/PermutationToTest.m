function test_suite = PermutationToTest()
    disp(['Setting up tests in ', mfilename()]);
    try
        test_functions = localfunctions();
    catch
    end
    initTestSuite;
end

function test_findPermutationsTo
    n = 10;
    d = 10;
    for i = 1:n
        S = replab.S(d);
        G = S.randomProperSubgroup(2);
        g = G.sample;
        s = randi(3, 1, d);
        t = s(g);
        P = G.findPermutationsTo(s, t);
        p = P.sample;
        assertEqual(s, t(g));
    end
end
