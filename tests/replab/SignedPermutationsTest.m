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
        G = replab.SignedSymmetricGroup(n);
        test_suite = replab.SignedPermutationGroupLaws(G).addTestCases(test_suite);
    end
end
