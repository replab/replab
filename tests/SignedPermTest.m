function test_suite = SignedPermTest()
    try
        test_functions = localfunctions();
    catch
    end
    initTestSuite;
    n = 10;
    G = replab.SignedPermutations(n);
    A = G.naturalAction;
    M = G.permutationIsomorphism;
    test_suite = G.lawsAddTestCases(test_suite);
    test_suite = A.lawsAddTestCases(test_suite, 'name', 'natural action');
    test_suite = M.lawsAddTestCases(test_suite, 'name', 'permutation isomorphism');
end
