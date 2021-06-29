function test_suite = Issue14Test()
     disp(['Setting up tests in ', mfilename()]);
    try
        test_functions = localfunctions();
    catch
    end
    initTestSuite;

    if ReplabTestParameters.onlyFastTests
        return;
    end

    group = replab.SignedPermutationGroup.of([1 -3 2 4]);
    I = group.naturalRep.decomposition;
    test_suite = I.laws.addTestCases(test_suite);
end
