function test_suite = WreathProductGroupTest()
    disp(['Setting up tests in ', mfilename()]);
    try
        test_functions = localfunctions();
    catch
    end
    initTestSuite;

    S3 = replab.S(2);
    W = S3.wreathProduct(S3);
    test_suite = W.laws.addTestCases(test_suite);
end

function test_wreath_subgroup
    S2 = replab.S(2);
    W = S2.wreathProduct(S2);
    g = {[2 1] {[1 2] [2 1]}};
    assertEqual(W.subgroup({g}).order, vpi(4));
end
