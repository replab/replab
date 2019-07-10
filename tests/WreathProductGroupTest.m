function test_suite = WreathProductGroupTest()
    try
        test_functions = localfunctions();
    catch
    end
    initTestSuite;
    S3 = replab.Permutations(2);
    W = replab.WreathProductGroup(S3, S3);
    test_suite = replab.FiniteGroupLaws(W).addTestCases(test_suite);
end
