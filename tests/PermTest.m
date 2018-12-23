function test_suite = PermTest()
    try
        test_functions = localfunctions();
    catch
    end
    initTestSuite;
    n = 10;
    G = replab.Permutations(n);
    A = G.naturalAction;
    test_suite = G.lawsAddTestCases(test_suite);
    test_suite = A.lawsAddTestCases(test_suite, 'name', 'natural action');
end
