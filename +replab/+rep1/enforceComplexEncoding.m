function W = enforceComplexEncoding(rep)
% Finds change of basis that reveals the complex division algebra
%
% Returns W such that W' * rep.image(g) * W is encoded using the
% quaternion encoding of replab.DivisionAlgebra
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
    W1 = w1;
    W2 = w2;
    while size(W1, 2)*2 < d
        A = rep.sample;
        if rank([W1 W2 A*w1]) > size(W1, 2)*2
            W1 = [W1 A*w1];
            W2 = [W2 A*w2];
        end
    end
    [Q R] = qr(W1, 0);
    W1 = W1 / R;
    W2 = W2 / R;
    W = reshape([W1; W2], [d d]);
end
