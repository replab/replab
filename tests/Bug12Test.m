function test_suite = Bug12Test()
    disp(['Setting up tests in ', mfilename()]);
    try
        test_functions = localfunctions();
    catch
    end
    initTestSuite;

    if ReplabTestParameters.onlyFastTests
        return;
    end

    S9 = replab.S(9);
    group = S9.subgroup({[7 4 1 9 5 2 6 8 3] [7 3 4 2 5 6 9 8 1]});
    I = group.definingRep.decomposition;
    test_suite = replab.IrreducibleLaws(I).addTestCases(test_suite);
end
