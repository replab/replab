function test_suite = PermGrpTest()
    try
        test_functions = localfunctions();
    catch
    end
    initTestSuite;
end

function test_factorization()
    n = 5;
    G = replab.PermGrp.symmetric(n);
    for i = 1:100
        p = randperm(n);
        w = G.factorization(p);
        p1 = G.evaluateWord(w);
        assertEqual(p1, p1);
    end
end

function test_randomBag()
    n = 5;
    G = replab.PermGrp.symmetric(n);
    for i = 1:100
        p = G.randomElement;
        assertTrue(G.contains(p));
    end
end
