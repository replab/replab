function test_suite = involutionTest()
    disp(['Setting up tests in ', mfilename()]);
    try
        test_functions = localfunctions();
    catch
    end
    initTestSuite;
end
function test_random
    A = randn(7,7) + 1i*randn(7,7);
    A = conj(A)*inv(A);
    [V, D, n] = replab.numerical.antilinear.decomposeInvolution(A);
    tol = 1e-13;
    assert(norm(A*V - V*D) < tol);
end
