function test_suite = AutomorphismGroup()
    disp(['Setting up tests in ', mfilename()]);
    try
        test_functions = localfunctions();
    catch
    end
    initTestSuite;

    S = replab.S(5);
    A = replab.AutomorphismGroup(S);
    test_suite = A.laws.addTestCases(test_suite);
end
