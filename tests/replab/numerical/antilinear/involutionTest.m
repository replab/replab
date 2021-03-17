function test_suite = involutionTest()
    disp(['Setting up tests in ', mfilename()]);
    try
        test_functions = localfunctions();
    catch
    end
    initTestSuite;
end
function test_random
    A = randn(25,25) + 1i*randn(25,25);
    A = conj(A)*inv(A);
    [V, D, n] = replab.numerical.antilinear.eigInvolution(A);
    tol = 1e-12;
    assert(norm(A*V - V*D) < tol);
end
