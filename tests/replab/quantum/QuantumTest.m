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
    if exist('syms') && ~replab.compat.isOctave
        %        G = replab.quantum.GeneralizedPauli(3);
        %test_suite = replab.FiniteGroupLaws(G).addTestCases(test_suite);
        %test_suite = replab.RepLaws(G.definingRep).addTestCases(test_suite);
        [G rep] = replab.quantum.clifford_qudit(2);
        test_suite = G.laws.addTestCases(test_suite);
        test_suite = rep.laws.addTestCases(test_suite);
    end

%     % The following test is currently commented as of issue #247
%     [G rep] = replab.quantum.clifford_qudit(3);
%     test_suite = replab.FiniteGroupLaws(G).addTestCases(test_suite);
%     test_suite = replab.RepLaws(rep).addTestCases(test_suite);
end
