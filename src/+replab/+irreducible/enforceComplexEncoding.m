function W = enforceComplexEncoding(rep, samples, sub)
% Finds change of basis that expresses the canonical basis of a complex division algebra
%
% See `replab.domain.ComplexTypeMatrices`
%
% Args:
%   rep (replab.Rep): Real representation being decomposed
%   samples (replab.irreducible.OnDemandSamples): Lazy evaluation of various samples for `rep`
%   sub (row cell array of replab.SubRep): Irreducible complex-type real subrepresentation of `rep`
%
% Returns
% -------
%   M: double matrix
%     Matrix `W` such that ``W * rep.image(g) * W'`` is in the canonical basis
    assert(isa(rep, 'replab.Rep'));
    assert(rep == sub.parent);
    assert(isequal(rep.field, 'R'));
    d = sub.dimension;
    S = sub.U*samples.commutantSample(2)*sub.U';
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
        g = rep.group.sample;
        x1 = sub.matrixRowAction(g, w1);
        if norm(x1 - W * (W' * x1)) > tol
            x2 = sub.matrixRowAction(g, w2);
            X = [x1 x2];
            X = X - W*(W'*X);
            t = trace(X'*X)/2;
            X = X/sqrt(t);
            W = [W X];
        end
    end
    W = W'; % returns adjoint because we piled up column vectors
end
