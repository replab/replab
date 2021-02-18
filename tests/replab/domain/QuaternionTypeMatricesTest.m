function test_suite = QuaternionTypeMatricesTest()
    disp(['Setting up tests in ', mfilename()]);
    try
        test_functions = localfunctions();
    catch
    end
    initTestSuite;
end

function test_back_and_forth
    A = randn(2,2); B = randn(2,2); C = randn(2,2); D = randn(2,2);
    X1 = replab.domain.QuaternionTypeMatrices.toMatrix(A, B, C, D, 'group');
    X2 = replab.domain.QuaternionTypeMatrices.toMatrix(A, B, C, D, 'commutant');
    [A1, B1, C1, D1] = replab.domain.QuaternionTypeMatrices.fromMatrix(X1, 'group');
    [A2, B2, C2, D2] = replab.domain.QuaternionTypeMatrices.fromMatrix(X2, 'commutant');
    assert(all(all(A == A1)) && all(all(A == A2)));
    assert(all(all(B == B1)) && all(all(B == B2)));
    assert(all(all(C == C1)) && all(all(C == C2)));
    assert(all(all(D == D1)) && all(all(D == D2)));
end
