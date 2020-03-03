function test_suite = DirectProductGroupTest()
    disp(['Setting up tests in ', mfilename()]);
    try
        test_functions = localfunctions();
    catch
    end
    initTestSuite;

    S4 = replab.Permutations(4);
    S4xS4 = S4.directProduct(S4);
    test_suite = replab.FiniteGroupLaws(S4xS4).addTestCases(test_suite);
end
