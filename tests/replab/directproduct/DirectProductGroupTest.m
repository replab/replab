function test_suite = DirectProductGroupTest()
    disp(['Setting up tests in ', mfilename()]);
    try
        test_functions = localfunctions();
    catch
    end
    initTestSuite;

    S4 = replab.S(4);
    S4xS4 = S4.directProduct(S4);
    test_suite = S4xS4.laws.addTestCases(test_suite);
end
