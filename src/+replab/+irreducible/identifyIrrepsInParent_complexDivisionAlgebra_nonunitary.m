function irreps = identifyIrrepsInParent_complexDivisionAlgebra_nonunitary(sub, sample)
% Identifies the irreducible representation(s) present in a real subrepresentation encoding a pair of conjugate irreps
%
% Args:
%   sub (`+replab.SubRep`): Subrepresentation with `+replab.Rep.divisionAlgebraName` set to ``'complex'``
%   sample (double(\*,\*)): Sample of ``sub.parent.commutant``
%
% Returns:
%   cell(1,\*) of `+replab.SubRep`: Real irreps with `+replab.Rep.frobeniusSchurIndicator` computed and `+replab.Rep.divisionAlgebraName` set
    assert(sub.overR);
    assert(strcmp(sub.divisionAlgebraName, 'complex'));
    d = sub.dimension;
    I = sub.injection('double/sparse');
    P = sub.projection('double/sparse');
    d = sub.dimension;
    A = I(:, 1:2:d);
    B = I(:, 2:2:d);
    C = P(1:2:d, :);
    D = P(2:2:d, :);
    U = [A+1i*B A-1i*B];
    V = [C-1i*D;C+1i*D];
    S = sample;
    X = V*S*U;
    E = diag(X);
    tol = 1e-10;
    if any(abs(E(1:d/2) - conj(E(d/2+1:end))) > tol)
        irreps = {};
        return
    end
    J = X(1:d/2, d/2+1:d);
    cJ = X(d/2+1:d, 1:d/2);
    cJJ = conj(J)*J;
    lambda = trace(cJJ)/(d/2);
    err = cJJ - lambda*eye(d/2);
    tol = replab.globals.doubleEigTol;
    if abs(lambda) < tol
        % complex-type representation
        sub.cache('frobeniusSchurIndicator', 0, '==');
        sub.cache('isIrreducible', true, '==');
        irreps = {sub};
    elseif lambda > 0
        [sub1, sub2] = replab.irreducible.regularizeRealPair(sub, sample);
        irreps = {sub1, sub2};
    else % lambda < 0
         % quaternion-type representation
        assert(mod(d, 4) == 0);
        irreps = {replab.irreducible.regularizeQuaternionic(sub, sample)};
    end
end
