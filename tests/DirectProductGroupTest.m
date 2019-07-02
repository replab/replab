function test_suite = DirectProductGroupTest()
    try
        test_functions = localfunctions();
    catch
    end
    initTestSuite;
    S4 = replab.Permutations(4);
    S4xS4 = replab.DirectProductGroup({S4 S4});
    test_suite = replab.FiniteGroupLaws(S4xS4).addTestCases(test_suite);
end
