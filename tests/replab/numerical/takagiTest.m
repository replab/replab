function test_suite = takagiTest()
    disp(['Setting up tests in ', mfilename()]);
    try
        test_functions = localfunctions();
    catch
    end
    initTestSuite;
end
function test_random
% We test the three algorithms for the Takagi decomposition
    n = 20; % size of the matrix
    % sample a complex symmetric matrix
    X = randn(n,n) + 1i*randn(n,n);
    J = X.'+X; % J is symmetric
    [U1,D1] = replab.numerical.takagi(J, 'svd');
    [U2,D2] = replab.numerical.takagi(J, 'hybrid');
    [U3,D3] = replab.numerical.takagi(J, 'jacobi');
    tol = 1e-10;
    assert(norm(U1.'*D1*U1 - J, 'fro') < tol);
    assert(norm(U2.'*D2*U2 - J, 'fro') < tol);
    assert(norm(U3.'*D3*U3 - J, 'fro') < tol);
end
