function res = enforceComplexEncoding(rep, context)
% Finds the similar representation that expresses the canonical basis of representation on a complex division algebra
%
%
% See `+replab.+domain.ComplexTypeMatrices`
%
% Args:
%   rep (`+replab.Rep`): Real representation with `rep.frobeniusSchurIndicator == 0`
%   context (`+replab.Context`): Sampling context
%
% Returns:
%   `+replab.SimilarRep`: Similar representation that has the proper encoding
    assert(isa(rep, 'replab.Rep'));
    assert(rep.overR);
    assert(isequal(rep.frobeniusSchurIndicator, 0));
    if isequal(rep.isDivisionAlgebraCanonical, true)
        res = replab.SimilarRep.identical(rep);
        return
    end
    if ~isequal(rep.isUnitary, true)
        res = replab.irreducible.enforceComplexEncoding(rep.unitarize, context);
        res = replab.rep.collapse(res);
        return
    end
    d = rep.dimension;
    S = rep.commutant.sampleInContext(context, 1);
    A = (S + S') + 1i * (S - S');
    v1 = replab.domain.Vectors('C', d).sample;
    v1 = v1/norm(v1);
    w1p = A*v1;
    alpha1 = real(w1p'*v1);
    w1 = w1p - alpha1*v1;
    beta2 = norm(w1);
    v2 = w1/beta2;
    w2p = A*v2;
    alpha2 = real(w2p'*v2);
    w2 = w2p - alpha2*v2 - beta2*v1;
    assert(norm(w2) < replab.Parameters.doubleEigTol);
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
    tol = replab.Parameters.doubleEigTol;
    while size(W, 2) < d
        g = rep.group.sample;
        x1 = rep.matrixRowAction(g, w1);
        if norm(x1 - W * (W' * x1)) > tol
            x2 = rep.matrixRowAction(g, w2);
            X = [x1 x2];
            X = X - W*(W'*X);
            t = trace(X'*X)/2;
            X = X/sqrt(t);
            W = [W X];
        end
    end
    res = rep.similarRep(W', W);
    res.isDivisionAlgebraCanonical = true;
end
