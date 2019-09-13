function W = enforceQuaternionEncoding(rep)
% Finds change of basis that reveals the quaternion division algebra
%
% Returns W such that W * rep.image(g) * W' is encoded using the
% quaternion encoding of replab.DivisionAlgebra
    assert(isa(rep, 'replab.Rep'));
    assert(isequal(rep.field, 'R'));
    d = rep.dimension;
    A1 = rep.commutant.sampleSelfAdjoint;
    Ai = rep.commutant.sampleSkewAdjoint;
    Aj = rep.commutant.sampleSkewAdjoint;
    Ak = rep.commutant.sampleSkewAdjoint;
    A = replab.quaternion.H(A1, Ai, Aj, Ak);
    v1 = replab.quaternion.Vectors(d).sample;
    v1 = v1/norm(v1);
    w1p = A*v1;
    alpha1 = real(w1p'*v1);
    w1 = w1p - alpha1*v1;
    beta2 = norm(w1);
    v2 = w1/beta2;
    w2p = A*v2;
    alpha2 = real(w2p'*v2);
    w2 = w2p - alpha2*v2 - beta2*v1;
    beta3 = norm(w2);
    v3 = w2/beta3;
    w3p = A*v3;
    alpha3 = real(w3p'*v3);
    w3 = w3p - alpha3*v3 - beta3*v2;
    beta4 = norm(w3);
    v4 = w3/beta4;
    w4p = A*v4;
    alpha4 = real(w4p'*v4);
    w4 = w4p - alpha4*v4 - beta4*v3;
    assert(norm(w4) < replab.Settings.doubleEigTol);
    T = [alpha1 beta2 0     0
         beta2 alpha2 beta3 0
         0     beta3  alpha3 beta4
         0      0     beta4  alpha4];
    V = [v1 v2 v3 v4];
    [U D] = eig(T);
    W = V*U;
    W1 = W(:,1);
    w1 = W1.r;
    w2 = W1.i;
    w3 = W1.j;
    w4 = W1.k;
    w1 = w1/norm(w1);
    w2 = w2/norm(w2);
    w3 = w3/norm(w3);
    w4 = w4/norm(w4);
    W = [w1 w2 w4 w3]; % switch to force basis (TODO: prove)
    % W is orthonormal
    tol = replab.Settings.doubleEigTol;
    while size(W, 2) < d
        A = rep.sample;
        x1 = A * w1;
        if norm(x1 - W * (W' * x1)) > tol
            x2 = A*w2;
            x3 = A*w3;
            x4 = A*w4;
            X = [x1 x2 x4 x3];
            X = X - W*(W'*X);
            t = trace(X'*X)/4;
            X = X/sqrt(t);
            % TODO: remove
            assert(~replab.isNonZeroMatrix(X'*X - eye(4), tol));
            W = [W X];
        end
    end
    W = W'; % returns adjoint as per new convention
end
