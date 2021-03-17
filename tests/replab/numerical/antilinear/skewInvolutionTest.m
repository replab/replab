function test_suite = skewInvolutionTest()
    disp(['Setting up tests in ', mfilename()]);
    try
        test_functions = localfunctions();
    catch
    end
    initTestSuite;
end
function test_random_non_symmetric
    d = 10;
    A = randn(2*d,2*d) + 1i*randn(2*d,2*d);
    J = [zeros(d,d) -eye(d); eye(d) zeros(d,d)];
    J = A*J*conj(inv(A));
    [V, D, n] = replab.numerical.antilinear.eigSkewInvolution(J);
    tol = 1e-12;
    assert(norm(J*V - V*D) < tol);
end
function test_random_symmetric
    d = 10;
    G = replab.U(d*2);
    U = G.sample;
    J = [zeros(d,d) -eye(d); eye(d) zeros(d,d)];
    J = U*J*U.';
    J = (J-J.')/2;
    [V, D, n] = replab.numerical.antilinear.eigSkewInvolution(J);
    tol = 1e-12;
    assert(norm(J*V - V*D) < tol);
end
