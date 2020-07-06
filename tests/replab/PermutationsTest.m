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
        G = replab.S(n);
        test_suite = replab.PermutationGroupLaws(G).addTestCases(test_suite);
    end
    S5 = replab.S(5);
    signRep = S5.signRep;
    test_suite = replab.RepLaws(signRep).addTestCases(test_suite);
end

function test_sign
    sigma = [3 4 5 2 1];
    assertEqual(replab.Permutation.sign(sigma), -1);
end
