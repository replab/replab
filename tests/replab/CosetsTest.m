function test_suite = CosetsTest()
    disp(['Setting up tests in ', mfilename()]);
    try
        test_functions = localfunctions();
    catch
    end
    initTestSuite;

    D8 = replab.DihedralGroup(8);
    C8 = replab.CyclicGroup(8);
    leftCosets = D8/C8;
    rightCosets = C8\D8;
    test_suite = replab.PermutationGroupLeftCosetsLaws(leftCosets).addTestCases(test_suite);
    test_suite = replab.PermutationGroupRightCosetsLaws(rightCosets).addTestCases(test_suite);
end
