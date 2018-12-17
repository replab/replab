function test_suite = WordTest()
    try
        test_functions = localfunctions();
    catch
    end
    initTestSuite;
    nG = 5;
    G = replab.FreeGroup(nG);
    test_suite = G.lawsAddTestCases(test_suite);
end
