function test_suite = LexMinTest()
    disp(['Setting up tests in ', mfilename()]);
    try
        test_functions = localfunctions();
    catch
    end
    initTestSuite;
end
function test_lexmin
    for i = 1:10
        n = 10;
        G = replab.S(n);
        G = G.randomProperSubgroup(3);
        s1 = randi(3, 1, n);
        g = G.sample;
        s2 = s1(g);
        [smin1, P1] = G.vectorFindLexMinimal(s1);
        [smin2, P2] = G.vectorFindLexMinimal(s2);
        assert(all(smin1 == smin2));
        assert(all(smin1(P1.sample) == s1));
        assert(all(smin2(P2.sample) == s2));
    end
end
function test_symmetric_group
    for i = 1:10
        n = 10;
        G = replab.S(n);
        s = randi(3, 1, n);
        [smin, P] = G.vectorFindLexMinimal(s);
        assert(all(smin == sort(s)));
    end
end
