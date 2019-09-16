function test_suite = LawsTest()
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
        G = replab.Permutations(n);
        test_suite = replab.PermutationsLaws(G).addTestCases(test_suite);
    end
    for n = allNs
        G = replab.signed.Permutations(n);
        test_suite = replab.signed.PermutationsLaws(G).addTestCases(test_suite);
    end
end
