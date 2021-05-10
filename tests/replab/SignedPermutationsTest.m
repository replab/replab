function test_suite = SignedPermutationsTest()
    disp(['Setting up tests in ', mfilename()]);
    try
        test_functions = localfunctions();
    catch
    end
    initTestSuite;

    if ReplabTestParameters.onlyFastTests
        allNs = [10];
    else
        allNs = [0 1 2 10];
    end

    for n = allNs
        G = replab.SignedSymmetricGroup(n);
        test_suite = G.laws.addTestCases(test_suite);
    end
end


function test_niceFiniteGroup_implementations
    G = replab.SignedSymmetricGroup(3);
    rep = G.regularRep;
    assert(rep.dimension == 48);
    assert(~G.isSimple);
    G1 = G.normalClosure([-1 2 3]);
    assert(G1.order == 8);
end
