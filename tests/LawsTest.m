function test_suite = LawsTest()
    try
        test_functions = localfunctions();
    catch
    end
    initTestSuite;
    for n = [0 1 2 10]
        G = replab.Permutations(n);
        test_suite = replab.PermutationsLaws(G).addTestCases(test_suite);
    end
    for n = [0 1 2 10]
        G = replab.SignedPermutations(n);
        test_suite = replab.SignedPermutationsLaws(G).addTestCases(test_suite);
    end
    for n = [0 1 5]
        G = replab.FreeGroup(n);
        test_suite = replab.FinitelyGeneratedGroupLaws(G).addTestCases(test_suite);
    end
end
