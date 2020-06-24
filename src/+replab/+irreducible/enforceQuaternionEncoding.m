function res = enforceQuaternionEncoding(rep, context)
% Finds the similar representation that expresses the canonical basis of representation on a quaternion division algebra
%
% See `+replab.+domain.QuaternionTypeMatrices`
%
% Args:
%   rep (`+replab.Rep`): Real representation with ``rep.frobeniusSchurIndicator == -2``
%   context (`replab.Context`): Sampling context
%
% Returns:
%   `+replab.SimilarRep`: Similar representation that has the proper encoding
    assert(isa(rep, 'replab.Rep'));
    assert(rep.overR);
    assert(isequal(rep.frobeniusSchurIndicator, -2));
    if isequal(rep.isDivisionAlgebraCanonical, true)
        res = replab.SimilarRep.identical(rep);
        return
    end
    if ~isequal(rep.isUnitary, true)
        res = replab.irreducible.enforceQuaternionEncoding(rep.unitarize, context);
        res = replab.rep.collapse(res);
        return
    end
    d = rep.dimension;
    S1 = rep.commutant.sampleInContext(context, 1);
    S2 = rep.commutant.sampleInContext(context, 2);
    S3 = rep.commutant.sampleInContext(context, 3);
    A1 = S1 + S1';
    Ai = S1 - S1';
    Aj = S2 - S2';
    Ak = S3 - S3';
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
    assert(norm(w4) < replab.Parameters.doubleEigTol);
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
    tol = replab.Parameters.doubleEigTol;
    while size(W, 2) < d
        g = rep.group.sample;
        x1 = rep.matrixRowAction(g, w1);
        if norm(x1 - W * (W' * x1)) > tol
            x2 = rep.matrixRowAction(g, w2);
            x3 = rep.matrixRowAction(g, w3);
            x4 = rep.matrixRowAction(g, w4);
            X = [x1 x2 x4 x3];
            X = X - W*(W'*X);
            t = trace(X'*X)/4;
            X = X/sqrt(t);
            % TODO: remove
            assert(~replab.isNonZeroMatrix(X'*X - eye(4), tol));
            W = [W X];
        end
    end
    res = rep.similarRep(W', W);
    res.isDivisionAlgebraCanonical = true;
end
