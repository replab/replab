function test_suite = SignedPermutationsTest()
    disp(['Setting up tests in ', mfilename()]);
    try
        test_functions = localfunctions();
    catch
    end
    initTestSuite;

    if ReplabTestParameters.onlyFastTests
        allNs = [10];
    else
        allNs = [0 1 2 10];
    end

    for n = allNs
        G = replab.SignedPermutations(n);
        test_suite = replab.SignedPermutationsLaws(G).addTestCases(test_suite);
    end
end
