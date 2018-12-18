function test_suite = PermTest()
    try
        test_functions = localfunctions();
    catch
    end
    initTestSuite;
    n = 10;
    G = replab.Permutations(n);
    A = G.naturalAction;
    A1 = G.vectorAction('R15');
    A2 = G.selfAdjointMatrixAction('R15');
    R = G.naturalRepresentation;
    test_suite = G.lawsAddTestCases(test_suite);
    test_suite = A.lawsAddTestCases(test_suite, 'name', 'natural action');
    test_suite = A1.lawsAddTestCases(test_suite, 'name', 'vector action');
    test_suite = A2.lawsAddTestCases(test_suite, 'name', 'self-adjoint matrix action');
    test_suite = R.lawsAddTestCases(test_suite, 'name', 'natural representation');
end

% $$$ 
% $$$ function test_vectorAction()
% $$$ % test compatibility of action with group
% $$$     n = 100;
% $$$     for i = 1:100
% $$$         v = rand(n, 1);
% $$$         x = replab.Perm.random(n);
% $$$         y = replab.Perm.random(n);
% $$$         z = replab.Perm.compose(x, y);
% $$$         imgy = replab.Perm.vectorAction(y, v);
% $$$         imgxy = replab.Perm.vectorAction(x, imgy);
% $$$         imgz = replab.Perm.vectorAction(z, v);
% $$$         assertEqual(imgz, imgxy);
% $$$         Mz = replab.Perm.matrix(z);
% $$$         assertEqual(Mz * v, imgz);
% $$$     end
% $$$ end
% $$$ 
% $$$ function test_selfAdjointMatrixAction()
% $$$     n = 10;
% $$$     for i = 1:100
% $$$         M = rand(n, n);
% $$$         x = replab.Perm.random(n);
% $$$         y = replab.Perm.random(n);
% $$$         z = replab.Perm.compose(x, y);
% $$$         imgy = replab.Perm.selfAdjointMatrixAction(y, M);
% $$$         imgxy = replab.Perm.selfAdjointMatrixAction(x, imgy);
% $$$         imgz = replab.Perm.selfAdjointMatrixAction(z, M);
% $$$         assertEqual(imgz, imgxy);
% $$$         rhoz = replab.Perm.matrix(z);
% $$$         assertEqual(rhoz * M * rhoz', imgz);
% $$$     end
% $$$ end
