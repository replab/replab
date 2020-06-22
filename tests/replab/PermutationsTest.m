function test_suite = PermutationsTest()
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
end

function test_sign
    sigma = [3 4 5 2 1];
    assertEqual(replab.Permutations.sign(sigma), -1);
end
