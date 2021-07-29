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
function test_issue_486
    g1 = [7 6 5 9 4 2 1 10 3 8];
    g2 = [7 6 9 5 3 2 8 1 4 10];
    G = replab.PermutationGroup.of(g1, g2);
    g = [7 2 4 3 9 6 1 10 5 8];
    assert(G.contains(g));
    s1 = [1 3 1 1 2 2 3 2 3 2];
    s2 = [3 3 1 1 3 2 1 2 2 2];
    assert(all(s2 == s1(g)));
    r1 = G.vectorFindLexMinimal(s1);
    r2 = G.vectorFindLexMinimal(s2);
    assert(all(r1 == r2));
end
