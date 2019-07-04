function W = enforceQuaternionEncoding(rep)
% Finds change of basis that reveals the quaternion division algebra
%
% Returns W such that W' * rep.image(g) * W is encoded using the
% quaternion encoding of replab.DivisionAlgebrax
    assert(isequal(rep.field, 'R'));
    d = rep.dimension;
    A = rep.quaternionification.commutant.sampleSelfAdjoint;
    v1 = replab.domain.QuaternionVectors(d).sample;
    v1 = v1/norm(v1);
    w1p = A*v1;
    alpha1 = part(w1p'*v1, 1);
    w1 = w1p - alpha1*v1;
    beta2 = norm(w1);
    v2 = w1/beta2;
    w2p = A*v2;
    alpha2 = part(w2p'*v2, 1);
    w2 = w2p - alpha2*v2 - beta2*v1;
    beta3 = norm(w2);
    v3 = w2/beta3;
    w3p = A*v3;
    alpha3 = part(w3p'*v3, 1);
    w3 = w3p - alpha3*v3 - beta3*v2;
    beta4 = norm(w3);
    v4 = w3/beta4;
    w4p = A*v4;
    alpha4 = part(w4p'*v4, 1);
    w4 = w4p - alpha4*v4 - beta4*v3;
    assert(norm(w4) < replab.Settings.doubleEigTol);
    T = [alpha1 beta2 0     0
         beta2 alpha2 beta3 0
         0     beta3  alpha3 beta4
         0      0     beta4  alpha4];
    V = [v1 v2 v3 v4];
    [U D] = eig(T);
    W = V*U;
    w1 = part(W(:,1), 1);
    w2 = part(W(:,1), 2);
    w3 = part(W(:,1), 3);
    w4 = part(W(:,1), 4);
    w1 = w1/norm(w1);
    w2 = w2/norm(w2);
    w3 = w3/norm(w3);
    w4 = w4/norm(w4);
    W = [w1 w2 w3 w4];
    while size(W, 2) < d
        A = rep.sample;
        x1 = A * w1;
        if norm(x1 - W * (W \ x1)) > replab.Settings.doubleEigTol
            x2 = A*w2;
            x3 = A*w3;
            x4 = A*w4;
            W = [W x1 x2 x3 x4];
        end
    end
    [Q R] = qr(W, 0);
    W = W / R;
end
