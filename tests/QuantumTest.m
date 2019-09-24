function test_suite = QuantumTest()
    disp(['Setting up tests in ', mfilename()]);
    try
        test_functions = localfunctions();
    catch
    end
    initTestSuite;
    if ReplabTestParameters.onlyFastTests
        return
    end
    G = replab.quantum.GHZ(3,2);
    test_suite = replab.GroupLaws(G).addTestCases(test_suite);
    %test_suite = replab.RepLaws(G.naturalRep).addTestCases(test_suite); % slow test with Octave
    G = replab.quantum.GeneralizedPauli(3);
    test_suite = replab.FiniteGroupLaws(G).addTestCases(test_suite);
    test_suite = replab.RepLaws(G.naturalRep).addTestCases(test_suite);
    [G rep] = replab.quantum.clifford_qudit(2);
    test_suite = replab.FiniteGroupLaws(G).addTestCases(test_suite);
    test_suite = replab.RepLaws(rep).addTestCases(test_suite);
    [G rep] = replab.quantum.clifford_qudit(3);
    test_suite = replab.FiniteGroupLaws(G).addTestCases(test_suite);
    test_suite = replab.RepLaws(rep).addTestCases(test_suite);
end
