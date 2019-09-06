function W = enforceComplexEncoding(rep)
% Finds change of basis that reveals the complex division algebra
%
% Returns W such that W * rep.image(g) * W' is encoded using the
% complex encoding of replab.DivisionAlgebra
    assert(isequal(rep.field, 'R'));
    d = rep.dimension;
    A = rep.complexification.commutant.sampleSelfAdjoint;
    v1 = replab.domain.ComplexVectors(d).sample;
    v1 = v1/norm(v1);
    w1p = A*v1;
    alpha1 = real(w1p'*v1);
    w1 = w1p - alpha1*v1;
    beta2 = norm(w1);
    v2 = w1/beta2;
    w2p = A*v2;
    alpha2 = real(w2p'*v2);
    w2 = w2p - alpha2*v2 - beta2*v1;
    assert(norm(w2) < replab.Settings.doubleEigTol);
    T = [alpha1 beta2
         beta2 alpha2];
    V = [v1 v2];
    [U D] = eig(T);
    W = V*U;
    w1 = real(W(:,1));
    w2 = imag(W(:,1));
    w1 = w1/norm(w1);
    w2 = w2/norm(w2);
    W = [w1 w2];
    tol = replab.Settings.doubleEigTol;
    while size(W, 2) < d
        A = rep.sample;
        x1 = A * w1;
        if norm(x1 - W * (W' * x1)) > tol
            x2 = A * w2;
            X = [x1 x2];
            X = X - W*(W'*X);
            t = trace(X'*X)/2;
            X = X/sqrt(t);
            W = [W X];
        end
    end
    W = W'; % returns adjoint
end
