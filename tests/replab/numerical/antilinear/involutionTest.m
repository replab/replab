function test_suite = involutionTest()
    disp(['Setting up tests in ', mfilename()]);
    try
        test_functions = localfunctions();
    catch
    end
    initTestSuite;
end
function test_random_non_symmetric
    d = 25;
    A = randn(d,d) + 1i*randn(d,d);
    A = conj(A)*inv(A);
    [V, D, n] = replab.numerical.antilinear.eigInvolution(A);
    tol = 1e-11;
    assert(norm(A*V - V*D) < tol);
end
function test_random_symmetric
    d = 25;
    U = replab.U(d);
    S = U.sample;
    A = S*S.';
    A = (A+A.')/2;
    [V, D, n] = replab.numerical.antilinear.eigInvolution(A);
    tol = 1e-11;
    assert(norm(A*V - V*D) < tol);
end
