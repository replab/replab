function f = frobeniusSchurIndicator_complex_irreducible(rep)
% Computes the Frobenius-Schur indicator of an irreducible subrepresentation
%
% Args:
%   rep (`+replab.Rep`): Representation to compute the Frobenius-Schur indicator of
%
% Returns:
%   integer: Frobenius-Schur indicator
    assert(rep.isIrreducible);
    tol = replab.globals.doubleEigTol;
    A = rep.equivariantFrom(rep.conj, 'type', 'double');
    d = rep.dimension;
    J = A.sample;
    cJJ = conj(J)*J;
    lambda = trace(cJJ)/d;
    err = cJJ - lambda*eye(d);
    assert(norm(err, 'fro') <= tol);
    if abs(lambda) < tol
        f = 0; % inequivalent to the dual
    elseif lambda > 0
        f = 1; % real representation
    elseif lambda < 0
        f = -1; % quaternion representation
    end
end
