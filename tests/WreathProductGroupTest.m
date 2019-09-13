function test_suite = WreathProductGroupTest()
    try
        test_functions = localfunctions();
    catch
    end
    initTestSuite;
    S3 = replab.Permutations(2);
    W = S3.wreathProduct(S3);
    test_suite = replab.FiniteGroupLaws(W).addTestCases(test_suite);
end
