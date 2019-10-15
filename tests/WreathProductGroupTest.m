function test_suite = WreathProductGroupTest()
    disp(['Setting up tests in ', mfilename()]);
    try
        test_functions = localfunctions();
    catch
    end
    initTestSuite;
    
    S3 = replab.Permutations(2);
    W = S3.wreathProduct(S3);
    test_suite = replab.FiniteGroupLaws(W).addTestCases(test_suite);
end

function test_wreath_subgroup
    S2 = replab.Permutations(2);
    W = S2.wreathProduct(S2);
    g = {[2 1] {[1 2] [2 1]}};
    assertEqual(W.subgroup({g}).order, 4);
end
