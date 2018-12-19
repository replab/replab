function test_suite = SignedPermTest()
    try
        test_functions = localfunctions();
    catch
    end
    initTestSuite;
    n = 10;
    G = replab.SignedPermutations(n);
    A = G.naturalAction;
    %A1 = G.vectorAction('R15');
    %A2 = G.selfAdjointMatrixAction('R15');
    M = G.permutationIsomorphism;
    %    R = G.naturalRepresentation;
    test_suite = G.lawsAddTestCases(test_suite);
    test_suite = A.lawsAddTestCases(test_suite, 'name', 'natural action');
    %test_suite = A1.lawsAddTestCases(test_suite, 'name', 'vector action');
    %test_suite = A2.lawsAddTestCases(test_suite, 'name', 'self-adjoint matrix action');
    test_suite = M.lawsAddTestCases(test_suite, 'name', 'permutation isomorphism');
    %    test_suite = R.lawsAddTestCases(test_suite, 'name', 'natural representation');
end
